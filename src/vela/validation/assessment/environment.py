"""阶段三全酶实验布局与精修代表姿态的空间映射。"""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass
from pathlib import Path

import gemmi

from vela.config import AppConfig
from vela.core.provenance import (
    JsonValue,
    atomic_write_json,
    atomic_write_text,
    sha256_file,
    utc_now,
)
from vela.core.typed_data import object_list, object_mapping
from vela.preparation.receptors.cleaning.conformers import (
    resolve_alternate_conformation,
)
from vela.preparation.receptors.models import ReceptorDefinition
from vela.validation.bound_states.assets import PEPTIDE_CHAIN, RECEPTOR_CHAIN
from vela.validation.models import EnvironmentReference, ValidationError
from vela.validation.records import file_record, read_document, validate_record
from vela.validation.refinement.planning import (
    refinement_identity,
    verify_refinement_plan,
)

MODEL_RECEPTOR_CHAIN = "R"


@dataclass(frozen=True, slots=True)
class EnvironmentAsset:
    """一个经过装配和链身份核对的全酶参考。"""

    reference_id: str
    path: Path
    receptor_chain_id: str
    beta_chain_ids: tuple[str, ...]
    other_catalytic_chain_ids: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class EnvironmentMappingOutcome:
    """一次后验环境映射的报告规模。"""

    manifest_path: Path
    mapping_count: int


@dataclass(frozen=True, slots=True)
class _Representative:
    candidate_id: str
    receptor_id: str
    cluster_id: str
    target: str
    path: Path
    sha256: str


def _read_structure(path: Path, *, name: str) -> gemmi.Structure:
    if not path.is_file():
        raise ValidationError(f"{name} does not exist: {path}")
    try:
        structure = gemmi.read_structure(str(path))
    except RuntimeError as exc:
        raise ValidationError(f"invalid {name}: {path}") from exc
    if len(structure) != 1:
        raise ValidationError(f"{name} must contain one model: {path}")
    return structure


def _chain(structure: gemmi.Structure, chain_id: str, *, name: str) -> gemmi.Chain:
    matches = [chain for chain in structure[0] if chain.name == chain_id]
    if len(matches) != 1:
        raise ValidationError(f"{name} chain is missing or ambiguous: {chain_id}")
    return matches[0]


def _polymer_chain(
    *, config: AppConfig, source: gemmi.Chain, output_name: str, identity: str
) -> gemmi.Chain:
    output = gemmi.Chain(output_name)
    for source_residue in source:
        if source_residue.entity_type != gemmi.EntityType.Polymer:
            continue
        residue = source_residue.clone()
        resolve_alternate_conformation(
            receptor_id=identity,
            residue=residue,
            settings=config.preparation.altloc,
        )
        output.add_residue(residue)
    if not output:
        raise ValidationError(f"environment chain contains no polymer: {identity}")
    return output


def _assembly_structure(
    *,
    config: AppConfig,
    reference: EnvironmentReference,
    receptor: ReceptorDefinition,
) -> gemmi.Structure:
    raw_path = config.paths.data_dir / "receptors" / "raw" / f"{receptor.pdb_id}.cif"
    raw = _read_structure(raw_path, name="raw full-enzyme structure")
    assembly = next(
        (item for item in raw.assemblies if item.name == reference.assembly_id), None
    )
    if assembly is None:
        raise ValidationError(
            f"configured biological assembly is missing: {reference.reference_id}"
        )
    expected = {
        receptor.author_chain_id,
        *reference.beta_author_chain_ids,
        *reference.other_catalytic_author_chain_ids,
    }
    declared = {
        chain_id
        for generator in assembly.generators
        for chain_id in (*generator.chains, *generator.subchains)
    }
    if not expected.issubset(declared):
        raise ValidationError(
            f"environment chains are not all declared in assembly {reference.assembly_id}"
        )
    assembled = raw.clone()
    assembled.transform_to_assembly(
        reference.assembly_id,
        gemmi.HowToNameCopiedChain.Short,
    )
    output = gemmi.Structure()
    output.name = reference.reference_id
    output.add_model(gemmi.Model(1))
    for chain_id in (
        receptor.author_chain_id,
        *reference.beta_author_chain_ids,
        *reference.other_catalytic_author_chain_ids,
    ):
        output[0].add_chain(
            _polymer_chain(
                config=config,
                source=_chain(assembled, chain_id, name=reference.reference_id),
                output_name=chain_id,
                identity=f"{reference.reference_id}_{chain_id}",
            )
        )
    output.setup_entities()
    return output


def prepare_environment_references(config: AppConfig) -> tuple[EnvironmentAsset, ...]:
    """从只读原始 mmCIF 制备配置声明的完整全酶装配参考。"""
    receptor_by_id = {item.receptor_id: item for item in config.receptors}
    root = config.paths.data_dir / "validation" / "environments"
    assets: list[EnvironmentAsset] = []
    entries: list[dict[str, JsonValue]] = []
    for reference in config.validation.environment_references:
        receptor = receptor_by_id.get(reference.receptor_id)
        if receptor is None:
            raise ValidationError(
                f"environment receptor is not registered: {reference.receptor_id}"
            )
        structure = _assembly_structure(
            config=config, reference=reference, receptor=receptor
        )
        output_dir = root / reference.reference_id
        output_path = output_dir / "assembly_reference.cif"
        atomic_write_text(output_path, structure.make_mmcif_document().as_string())
        asset = EnvironmentAsset(
            reference_id=reference.reference_id,
            path=output_path,
            receptor_chain_id=receptor.author_chain_id,
            beta_chain_ids=reference.beta_author_chain_ids,
            other_catalytic_chain_ids=reference.other_catalytic_author_chain_ids,
        )
        assets.append(asset)
        entries.append(
            {
                "reference_id": reference.reference_id,
                "pdb_id": receptor.pdb_id,
                "assembly_id": reference.assembly_id,
                "receptor_id": receptor.receptor_id,
                "receptor_author_chain_id": receptor.author_chain_id,
                "beta_author_chain_ids": list(reference.beta_author_chain_ids),
                "other_catalytic_author_chain_ids": list(
                    reference.other_catalytic_author_chain_ids
                ),
                "evaluation_targets": list(reference.evaluation_targets),
                "construct_note": reference.construct_note,
                "output": file_record(output_path, root=root),
            }
        )
    atomic_write_json(
        root / "preparation_manifest.json",
        {
            "schema": "vela.environment-preparation-manifest/1",
            "generated_at": utc_now(),
            "evidence_category": "full_enzyme_environment_reference",
            "known_site_information_used": True,
            "entries": entries,
        },
    )
    return tuple(assets)


def _residue_key(residue: gemmi.Residue) -> tuple[int, str]:
    number = residue.seqid.num
    if number is None:
        raise ValidationError("environment alignment residue lacks a sequence number")
    return number, residue.seqid.icode.strip()


def _ca_residues(chain: gemmi.Chain) -> tuple[gemmi.Residue, ...]:
    return tuple(residue for residue in chain if _named_atom(residue, "CA") is not None)


def _named_atom(residue: gemmi.Residue, name: str) -> gemmi.Atom | None:
    for atom in residue:
        if atom.name == name:
            return atom
    return None


def _alignment(
    *, fixed: gemmi.Chain, moving: gemmi.Chain
) -> tuple[gemmi.Transform, float, int]:
    moving_by_key = {_residue_key(residue): residue for residue in _ca_residues(moving)}
    fixed_positions: list[gemmi.Position] = []
    moving_positions: list[gemmi.Position] = []
    for fixed_residue in _ca_residues(fixed):
        moving_residue = moving_by_key.get(_residue_key(fixed_residue))
        if moving_residue is None or moving_residue.name != fixed_residue.name:
            continue
        fixed_atom = _named_atom(fixed_residue, "CA")
        moving_atom = _named_atom(moving_residue, "CA")
        if fixed_atom is None or moving_atom is None:
            continue
        fixed_positions.append(fixed_atom.pos)
        moving_positions.append(moving_atom.pos)
    if len(fixed_positions) < 50:
        raise ValidationError(
            "fewer than 50 identity-matched CA atoms are available for environment alignment"
        )
    result = gemmi.superpose_positions(fixed_positions, moving_positions)
    return result.transform, result.rmsd, len(fixed_positions)


def _transformed_chain(
    source: gemmi.Chain, *, name: str, transform: gemmi.Transform
) -> gemmi.Chain:
    chain = source.clone()
    chain.name = name
    for residue in chain:
        for atom in residue:
            atom.pos = gemmi.Position(transform.apply(atom.pos))
    return chain


def _heavy_atoms(chains: tuple[gemmi.Chain, ...]) -> tuple[gemmi.Atom, ...]:
    return tuple(
        atom
        for chain in chains
        for residue in chain
        for atom in residue
        if atom.element.name != "H"
    )


def _spatial_metrics(
    *, peptide: gemmi.Chain, context: tuple[gemmi.Chain, ...], contact_A: float
) -> tuple[int, float]:
    peptide_atoms = _heavy_atoms((peptide,))
    context_atoms = _heavy_atoms(context)
    if not peptide_atoms or not context_atoms:
        raise ValidationError("environment mapping lacks heavy atoms")
    pairs = 0
    minimum = float("inf")
    for peptide_atom in peptide_atoms:
        for context_atom in context_atoms:
            distance = peptide_atom.pos.dist(context_atom.pos)
            minimum = min(minimum, distance)
            if distance <= contact_A:
                pairs += 1
    return pairs, minimum


def _representatives(
    *, run_dir: Path, plan: dict[str, object]
) -> tuple[_Representative, ...]:
    analysis_dir = run_dir / "refinement_analysis"
    manifest = read_document(
        analysis_dir / "analysis_manifest.json", name="refinement analysis manifest"
    )
    evidence_category, known_site_information_used = refinement_identity(plan)
    if (
        manifest.get("schema") != "vela.validation-refinement-analysis-manifest/3"
        or manifest.get("status") != "completed"
        or manifest.get("evidence_category") != evidence_category
        or manifest.get("known_site_information_used")
        is not known_site_information_used
    ):
        raise ValidationError("refinement analysis identity is invalid")
    clusters_path, _ = validate_record(
        root=analysis_dir,
        raw=manifest.get("refined_clusters"),
        name="refined clusters",
    )
    target_by_identity: dict[tuple[str, str], str] = {}
    try:
        tasks = object_list(plan.get("tasks"), name="refinement plan tasks")
    except TypeError as exc:
        raise ValidationError("refinement plan tasks are invalid") from exc
    for raw in tasks:
        try:
            task = object_mapping(raw, name="refinement plan task")
        except TypeError as exc:
            raise ValidationError("refinement plan task is invalid") from exc
        candidate = task.get("candidate_id")
        receptor = task.get("receptor_id")
        target = task.get("target")
        if (
            not isinstance(candidate, str)
            or not isinstance(receptor, str)
            or not isinstance(target, str)
        ):
            raise ValidationError("refinement plan task identity is invalid")
        target_by_identity[(candidate, receptor)] = target
    with clusters_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "cluster_id",
            "candidate_id",
            "receptor_id",
            "representative_path",
            "supported",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValidationError("refined cluster table columns are invalid")
        rows = tuple(dict(row) for row in reader)
    result: list[_Representative] = []
    for row in rows:
        if row["supported"] != "true":
            continue
        relative = Path(row["representative_path"])
        path = (run_dir / relative).resolve()
        try:
            path.relative_to(run_dir.resolve())
        except ValueError as exc:
            raise ValidationError("representative path escapes refinement run") from exc
        identity = (row["candidate_id"], row["receptor_id"])
        target = target_by_identity.get(identity)
        if target is None or not path.is_file():
            raise ValidationError("supported representative identity is invalid")
        result.append(
            _Representative(
                candidate_id=identity[0],
                receptor_id=identity[1],
                cluster_id=row["cluster_id"],
                target=target,
                path=path,
                sha256=sha256_file(path),
            )
        )
    if not result:
        raise ValidationError("refinement analysis has no supported representatives")
    return tuple(result)


def _tsv(rows: list[dict[str, str]]) -> str:
    fields: list[str] = [
        "mapping_id",
        "candidate_id",
        "receptor_id",
        "cluster_id",
        "target",
        "reference_id",
        "layout_kind",
        "matched_ca_atoms",
        "alignment_rmsd_A",
        "beta_contact_pairs",
        "minimum_beta_distance_A",
        "other_catalytic_contact_pairs",
        "minimum_other_catalytic_distance_A",
        "model_path",
        "model_sha256",
    ]
    output = io.StringIO()
    writer = csv.DictWriter(
        output, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def map_refinement_environment(
    *, config: AppConfig, refinement_run_dir: Path
) -> EnvironmentMappingOutcome:
    """将受支持精修代表姿态叠合到全酶布局并报告原始空间量。"""
    output_dir = refinement_run_dir / "environment_mapping"
    if output_dir.exists():
        raise ValidationError(f"environment mapping already exists: {output_dir}")
    plan, _ = verify_refinement_plan(config=config, run_dir=refinement_run_dir)
    source_evidence, source_known = refinement_identity(plan)
    representatives = _representatives(run_dir=refinement_run_dir, plan=plan)
    receptor_by_id = {item.receptor_id: item for item in config.receptors}
    root = config.paths.data_dir / "validation" / "environments"
    rows: list[dict[str, str]] = []
    index = 1
    for representative in representatives:
        candidate = _read_structure(representative.path, name="refined representative")
        moving_receptor = _chain(candidate, RECEPTOR_CHAIN, name="candidate receptor")
        moving_peptide = _chain(candidate, PEPTIDE_CHAIN, name="candidate peptide")
        for reference in config.validation.environment_references:
            if representative.target not in reference.evaluation_targets:
                continue
            receptor_definition = receptor_by_id[reference.receptor_id]
            reference_path = root / reference.reference_id / "assembly_reference.cif"
            reference_structure = _read_structure(
                reference_path, name="full-enzyme reference"
            )
            fixed_receptor = _chain(
                reference_structure,
                receptor_definition.author_chain_id,
                name="full-enzyme receptor",
            )
            transform, rmsd, matched = _alignment(
                fixed=fixed_receptor, moving=moving_receptor
            )
            model_receptor = _transformed_chain(
                moving_receptor, name=MODEL_RECEPTOR_CHAIN, transform=transform
            )
            peptide = _transformed_chain(
                moving_peptide, name=PEPTIDE_CHAIN, transform=transform
            )
            beta = tuple(
                _chain(reference_structure, chain_id, name="CK2beta context").clone()
                for chain_id in reference.beta_author_chain_ids
            )
            other = tuple(
                _chain(reference_structure, chain_id, name="catalytic context").clone()
                for chain_id in reference.other_catalytic_author_chain_ids
            )
            beta_pairs, beta_minimum = _spatial_metrics(
                peptide=peptide,
                context=beta,
                contact_A=config.validation.interface_contact_A,
            )
            other_pairs, other_minimum = _spatial_metrics(
                peptide=peptide,
                context=other,
                contact_A=config.validation.interface_contact_A,
            )
            mapping_id = f"environment_{index:04d}"
            index += 1
            model = gemmi.Structure()
            model.name = mapping_id
            model.add_model(gemmi.Model(1))
            model[0].add_chain(model_receptor)
            model[0].add_chain(peptide)
            for context_chain in (*beta, *other):
                model[0].add_chain(context_chain)
            model.setup_entities()
            model_path = output_dir / "models" / f"{mapping_id}.cif"
            atomic_write_text(model_path, model.make_mmcif_document().as_string())
            layout_kind = (
                "experimental_ck2alpha_layout_mapping"
                if representative.target == receptor_definition.target
                else "homology_ck2alpha_prime_layout_model"
            )
            rows.append(
                {
                    "mapping_id": mapping_id,
                    "candidate_id": representative.candidate_id,
                    "receptor_id": representative.receptor_id,
                    "cluster_id": representative.cluster_id,
                    "target": representative.target,
                    "reference_id": reference.reference_id,
                    "layout_kind": layout_kind,
                    "matched_ca_atoms": str(matched),
                    "alignment_rmsd_A": f"{rmsd:.6f}",
                    "beta_contact_pairs": str(beta_pairs),
                    "minimum_beta_distance_A": f"{beta_minimum:.6f}",
                    "other_catalytic_contact_pairs": str(other_pairs),
                    "minimum_other_catalytic_distance_A": f"{other_minimum:.6f}",
                    "model_path": model_path.relative_to(refinement_run_dir).as_posix(),
                    "model_sha256": sha256_file(model_path),
                }
            )
    if not rows:
        raise ValidationError("no representative matches an environment target")
    report_path = output_dir / "environment_mappings.tsv"
    atomic_write_text(report_path, _tsv(rows))
    analysis_manifest = (
        refinement_run_dir / "refinement_analysis" / "analysis_manifest.json"
    )
    reference_manifest = root / "preparation_manifest.json"
    manifest_path = output_dir / "mapping_manifest.json"
    atomic_write_json(
        manifest_path,
        {
            "schema": "vela.environment-mapping-manifest/1",
            "stage": "validation_environment_mapping",
            "status": "completed",
            "generated_at": utc_now(),
            "source_evidence_category": source_evidence,
            "source_known_site_information_used": source_known,
            "evidence_category": "post_hoc_full_enzyme_environment_mapping",
            "classification_applied": False,
            "parameters": {
                "contact_distance_A": config.validation.interface_contact_A,
                "steric_clash_threshold": "not_configured",
            },
            "inputs": {
                "refinement_analysis": file_record(
                    analysis_manifest, root=refinement_run_dir
                ),
                "environment_references": file_record(
                    reference_manifest, root=config.paths.data_dir
                ),
            },
            "mappings": {
                "path": report_path.name,
                "sha256": sha256_file(report_path),
                "count": len(rows),
            },
        },
    )
    return EnvironmentMappingOutcome(manifest_path, len(rows))
