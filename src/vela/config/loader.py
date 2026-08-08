"""按优先级组合默认值、项目 TOML 和环境覆盖。"""

from pathlib import Path

from vela.config.chemistry import parse_chemistry
from vela.config.design import parse_design
from vela.config.discovery import parse_discovery
from vela.config.models import AppConfig, ConfigError
from vela.config.preparation import parse_audit, parse_preparation
from vela.config.receptors import parse_receptors
from vela.config.sections import parse_download, parse_paths
from vela.config.validation import parse_validation
from vela.config.values import (
    apply_environment,
    assert_keys,
    default_document,
    merge,
    merge_disjoint,
    read_toml,
)
from vela.core.provenance import sha256_text

PROJECT_CONFIG_FILENAMES = (
    "common.toml",
    "discovery.toml",
    "validation.toml",
    "design.toml",
)


def _read_project_documents(
    config_dir: Path,
) -> tuple[dict[str, object], tuple[Path, ...], str]:
    """读取固定的四份项目参数文件并生成可追溯快照。"""
    if not config_dir.is_dir():
        raise ConfigError(f"config directory does not exist: {config_dir}")
    source_files = tuple(config_dir / name for name in PROJECT_CONFIG_FILENAMES)
    missing = [path.name for path in source_files if not path.is_file()]
    if missing:
        raise ConfigError("config directory is missing files: " + ", ".join(missing))

    external: dict[str, object] = {}
    snapshot_parts: list[str] = []
    for source_file in source_files:
        external = merge_disjoint(external, read_toml(source_file))
        content = source_file.read_text(encoding="utf-8")
        snapshot_parts.append(f"# --- {source_file.name} ---\n{content}")
        if not content.endswith("\n"):
            snapshot_parts.append("\n")
    return external, source_files, "".join(snapshot_parts)


def _validate_cross_section(config: AppConfig) -> None:
    """校验只有在完整项目配置中才能判断的跨 section 关系。"""
    if set(config.discovery.seeds) & set(config.discovery.qualification.seeds):
        raise ConfigError(
            "discovery production and qualification seeds must not overlap"
        )
    receptor_by_id = {item.receptor_id: item for item in config.receptors}
    for target in config.discovery.targets:
        for role, receptor_id in (
            ("reference", target.reference_receptor),
            ("pilot", target.pilot_receptor),
        ):
            receptor = receptor_by_id.get(receptor_id)
            if (
                receptor is None
                or receptor.target != target.target_id
                or "blind_discovery" not in receptor.roles
            ):
                raise ConfigError(
                    f"discovery target {target.target_id} {role} receptor must be a matching blind_discovery member"
                )
    state_by_id = {item.state_id: item for item in config.validation.bound_states}
    control_state = state_by_id.get(
        config.discovery.qualification.control_bound_state_id
    )
    if control_state is None or control_state.local_control_kind != (
        "standard_cyclic_peptide"
    ):
        raise ConfigError(
            "discovery qualification control must reference a standard cyclic peptide bound state"
        )
    native_receptor = receptor_by_id.get(control_state.receptor_id)
    control_receptors = tuple(
        receptor_by_id.get(receptor_id)
        for receptor_id in config.discovery.qualification.control_receptor_ids
    )
    benchmark_receptor = receptor_by_id.get(
        config.discovery.qualification.benchmark_receptor_id
    )
    configured_targets = {target.target_id for target in config.discovery.targets}
    if (
        config.discovery.qualification.control_target_id not in configured_targets
        or native_receptor is None
        or native_receptor.target != config.discovery.qualification.control_target_id
        or any(
            receptor is None
            or receptor.target != config.discovery.qualification.control_target_id
            or "blind_discovery" not in receptor.roles
            or "qualification_control" not in receptor.roles
            for receptor in control_receptors
        )
    ):
        raise ConfigError(
            "discovery qualification receptors must be matching blind-discovery controls"
        )
    if (
        benchmark_receptor is None
        or benchmark_receptor.target != config.discovery.qualification.control_target_id
        or "qualification_benchmark" not in benchmark_receptor.roles
    ):
        raise ConfigError(
            "discovery qualification benchmark receptor must match control_target_id"
        )
    if control_state.ligand_sequence is None or len(
        config.discovery.qualification.control_secondary_structure
    ) != len(control_state.ligand_sequence):
        raise ConfigError(
            "discovery qualification control secondary structure length must match its ligand"
        )
    topology_threshold = (
        config.discovery.cabsdock.max_reconstructable_disulfide_ca_distance_A
    )
    topology_candidates = (
        config.discovery.qualification.topology_calibration.candidate_ca_thresholds_A
    )
    if topology_threshold not in topology_candidates:
        raise ConfigError(
            "the CABS topology threshold must be one of the calibrated candidate thresholds"
        )
    for template in config.validation.guided_templates:
        if len(template.ligand_positions) != len(config.chemistry.sequence):
            raise ConfigError(
                f"guided template {template.template_id} must map every ligand residue"
            )
        state = state_by_id[template.bound_state_id]
        mapped_bonds = {
            tuple(
                sorted(
                    (
                        template.ligand_positions[bond.first - 1],
                        template.ligand_positions[bond.second - 1],
                    )
                )
            )
            for bond in config.chemistry.disulfide_bonds
        }
        state_bonds = {(bond.first, bond.second) for bond in state.disulfide_bonds}
        if not mapped_bonds or not mapped_bonds.issubset(state_bonds):
            raise ConfigError(
                f"guided template {template.template_id} does not preserve the ligand disulfide topology"
            )

    for reference in config.validation.environment_references:
        receptor = receptor_by_id.get(reference.receptor_id)
        if receptor is None:
            raise ConfigError(
                f"environment reference {reference.reference_id} has an unknown receptor"
            )
        if "full_enzyme_environment" not in receptor.roles:
            raise ConfigError(
                f"environment receptor {receptor.receptor_id} lacks its required role"
            )
        context_chains = (
            reference.beta_author_chain_ids + reference.other_catalytic_author_chain_ids
        )
        if receptor.author_chain_id in context_chains:
            raise ConfigError(
                f"environment reference {reference.reference_id} reuses its aligned receptor chain as context"
            )

    sequence_length = len(config.chemistry.sequence)
    disulfide_positions = {
        position
        for bond in config.chemistry.disulfide_bonds
        for position in (bond.first, bond.second)
    }
    if any(
        position > sequence_length or position in disulfide_positions
        for position in config.design.sequence.mutable_positions
    ):
        raise ConfigError(
            "design mutable_positions must be inside the ligand and exclude disulfide positions"
        )


def load_config(path: Path) -> AppConfig:
    """加载包默认值、四份外部 TOML 和环境变量并建立领域对象。"""
    source_dir = path.expanduser().resolve()
    external, source_files, snapshot_text = _read_project_documents(source_dir)
    document = apply_environment(merge(default_document(), external))
    assert_keys(
        document,
        allowed={
            "paths",
            "download",
            "audit",
            "preparation",
            "chemistry",
            "receptors",
            "discovery",
            "validation",
            "design",
        },
        required={
            "paths",
            "download",
            "audit",
            "preparation",
            "chemistry",
            "receptors",
            "discovery",
            "validation",
            "design",
        },
        path="root",
    )
    config = AppConfig(
        source_dir=source_dir,
        source_files=source_files,
        source_snapshot_text=snapshot_text,
        source_snapshot_sha256=sha256_text(snapshot_text),
        paths=parse_paths(document, config_dir=source_dir),
        download=parse_download(document),
        audit=parse_audit(document),
        preparation=parse_preparation(document),
        chemistry=parse_chemistry(document),
        receptors=parse_receptors(document),
        discovery=parse_discovery(document, config_dir=source_dir),
        validation=parse_validation(document, config_dir=source_dir),
        design=parse_design(document, config_dir=source_dir),
    )
    _validate_cross_section(config)
    return config
