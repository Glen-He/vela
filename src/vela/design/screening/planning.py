"""阶段四成对界面初筛的候选展开、resfile 和计划写入。"""

from __future__ import annotations

import csv
import shutil
from collections import defaultdict
from importlib.resources import files
from pathlib import Path

import gemmi

from vela.config import AppConfig
from vela.core.provenance import (
    JsonValue,
    atomic_write_json,
    atomic_write_text,
    sha256_file,
    sha256_text,
    utc_now,
    vela_software_identity,
)
from vela.core.run_identity import validate_run_id
from vela.core.typed_data import object_mapping
from vela.design.models import (
    DesignError,
    DesignTemplate,
    ScreenTask,
    SequenceCandidate,
)
from vela.design.readiness import assess_design_readiness
from vela.design.screening.inputs import selected_templates
from vela.design.screening.records import (
    SCREEN_PLAN_NAME,
    ScreenPlan,
    design_parameters,
    read_screen_plan,
)
from vela.design.sequence.library import (
    candidate_table,
    combination_library,
    systematic_single_library,
)
from vela.validation.records import (
    file_record,
    read_document,
    validate_record,
)
from vela.validation.rosetta import verify_rosetta_scripts_tool


def _polymer_residues(chain: gemmi.Chain) -> tuple[gemmi.Residue, ...]:
    return tuple(
        residue for residue in chain if any(atom.name == "CA" for atom in residue)
    )


def _residue_token(residue: gemmi.Residue, chain_id: str) -> tuple[int, str, str]:
    number = residue.seqid.num
    insertion = residue.seqid.icode.strip()
    if number is None:
        raise DesignError("design template contains a residue without a PDB number")
    return number, insertion, chain_id


def _token_text(token: tuple[int, str, str]) -> str:
    number, insertion, chain_id = token
    return f"{number}{insertion} {chain_id}"


def _neighbor_tokens(
    *, path: Path, mutation_positions: tuple[int, ...], distance_A: float
) -> tuple[tuple[int, str, str], ...]:
    try:
        structure = gemmi.read_structure(str(path))
    except RuntimeError as exc:
        raise DesignError(f"invalid design template: {path}") from exc
    if len(structure) != 1:
        raise DesignError("design template must contain one model")
    chains = {chain.name: chain for chain in structure[0]}
    if set(chains) != {"A", "P"}:
        raise DesignError("design template must contain only A/P chains")
    peptide = _polymer_residues(chains["P"])
    by_number = {residue.seqid.num: residue for residue in peptide}
    if len(by_number) != len(peptide):
        raise DesignError("ligand template residue numbering is duplicated")
    try:
        focus = tuple(by_number[position] for position in mutation_positions)
    except KeyError as exc:
        raise DesignError(
            "candidate mutation position is absent from its template"
        ) from exc
    focus_atoms = tuple(
        atom for residue in focus for atom in residue if atom.element.name != "H"
    )
    selected: set[tuple[int, str, str]] = set()
    for chain_id in ("A", "P"):
        for residue in _polymer_residues(chains[chain_id]):
            if any(
                atom.element.name != "H" and atom.pos.dist(focus_atom.pos) <= distance_A
                for atom in residue
                for focus_atom in focus_atoms
            ):
                selected.add(_residue_token(residue, chain_id))
    selected.update(_residue_token(residue, "P") for residue in focus)
    return tuple(sorted(selected, key=lambda item: (item[2] != "A", item[0], item[1])))


def _resfile_text(
    *,
    template: DesignTemplate,
    candidate: SequenceCandidate,
    reference_sequence: str,
    neighbor_distance_A: float,
    state: str,
) -> str:
    tokens = _neighbor_tokens(
        path=template.path,
        mutation_positions=candidate.mutation_positions,
        distance_A=neighbor_distance_A,
    )
    mutations = {
        position: (
            reference_sequence[position - 1]
            if state == "wt"
            else candidate.sequence[position - 1]
        )
        for position in candidate.mutation_positions
    }
    lines = ["NATRO\n", "start\n"]
    for token in tokens:
        number, insertion, chain_id = token
        residue = mutations.get(number) if chain_id == "P" and not insertion else None
        instruction = (
            "NATAA" if residue is None or state == "wt" else f"PIKAA {residue}"
        )
        lines.append(f"{_token_text(token)} {instruction}\n")
    return "".join(lines)


def _copy_templates(
    *, run_dir: Path, templates: tuple[DesignTemplate, ...]
) -> tuple[DesignTemplate, ...]:
    destination_dir = run_dir / "templates"
    destination_dir.mkdir(parents=True)
    copied: list[DesignTemplate] = []
    for template in templates:
        destination = destination_dir / f"{template.template_id}.pdb"
        shutil.copyfile(template.path, destination)
        digest = sha256_file(destination)
        if digest != template.sha256:
            raise DesignError("copied design template hash changed")
        copied.append(
            DesignTemplate(
                template_id=template.template_id,
                evidence_role=template.evidence_role,
                cluster_id=template.cluster_id,
                candidate_id=template.candidate_id,
                receptor_id=template.receptor_id,
                target=template.target,
                path=destination,
                sha256=digest,
                receptor_residue_count=template.receptor_residue_count,
                fixed_histidine_pose_indices=template.fixed_histidine_pose_indices,
            )
        )
    return tuple(copied)


def _write_protocol(run_dir: Path) -> Path:
    protocol = files("vela.resources").joinpath("design-interface-screen.xml")
    path = run_dir / "protocols" / "interface_screen.xml"
    atomic_write_text(path, protocol.read_text(encoding="utf-8"))
    return path


def _write_resfiles_and_tasks(
    *,
    config: AppConfig,
    run_dir: Path,
    candidates: tuple[SequenceCandidate, ...],
    templates: tuple[DesignTemplate, ...],
) -> tuple[ScreenTask, ...]:
    resfile_dir = run_dir / "resfiles"
    tasks: list[ScreenTask] = []
    task_index = 1
    pair_index = 1
    for candidate in candidates:
        for template in templates:
            paths: dict[str, tuple[Path, str]] = {}
            for state in ("wt", "mutant"):
                path = (
                    resfile_dir
                    / f"{candidate.candidate_id}__{template.template_id}__{state}.resfile"
                )
                atomic_write_text(
                    path,
                    _resfile_text(
                        template=template,
                        candidate=candidate,
                        reference_sequence=config.chemistry.sequence,
                        neighbor_distance_A=(config.design.screen.neighbor_distance_A),
                        state=state,
                    ),
                )
                paths[state] = path, sha256_file(path)
            for seed in config.design.seeds:
                pair_id = f"pair_{pair_index:07d}"
                pair_index += 1
                for state in ("wt", "mutant"):
                    path, digest = paths[state]
                    wt_context_id = (
                        "wt_"
                        + sha256_text(
                            "\0".join(
                                (
                                    template.template_id,
                                    str(seed),
                                    config.chemistry.chemistry_id,
                                    digest,
                                )
                            )
                        )[:16]
                        if state == "wt"
                        else None
                    )
                    tasks.append(
                        ScreenTask(
                            task_id=f"screen_{task_index:07d}",
                            pair_id=pair_id,
                            state=state,
                            candidate=candidate,
                            template=template,
                            seed=seed,
                            resfile_path=path,
                            resfile_sha256=digest,
                            wt_context_id=wt_context_id,
                        )
                    )
                    task_index += 1
    if not tasks:
        raise DesignError("design screen task matrix is empty")
    return tuple(tasks)


def candidate_record(candidate: SequenceCandidate) -> dict[str, JsonValue]:
    return {
        "candidate_id": candidate.candidate_id,
        "sequence": candidate.sequence,
        "mutation_string": candidate.mutation_string,
        "mutation_positions": list(candidate.mutation_positions),
        "mutation_count": candidate.mutation_count,
        "net_charge": candidate.net_charge,
        "design_round": candidate.design_round,
        "generation": candidate.generation,
        "parents": [
            {
                "candidate_id": item.candidate_id,
                "edit": item.edit,
                "edit_type": item.edit_type,
            }
            for item in candidate.parents
        ],
        "ancestor_candidate_ids": list(candidate.ancestor_candidate_ids),
        "proposal_source": candidate.proposal_source,
    }


def template_record(template: DesignTemplate, *, run_dir: Path) -> dict[str, JsonValue]:
    return {
        "template_id": template.template_id,
        "evidence_role": template.evidence_role,
        "cluster_id": template.cluster_id,
        "candidate_id": template.candidate_id,
        "receptor_id": template.receptor_id,
        "target": template.target,
        "structure": file_record(template.path, root=run_dir),
        "receptor_residue_count": template.receptor_residue_count,
        "fixed_histidine_pose_indices": list(template.fixed_histidine_pose_indices),
    }


def write_screen_plan(
    *,
    config: AppConfig,
    run_id: str,
    design_round: str,
    candidates: tuple[SequenceCandidate, ...],
    templates: tuple[DesignTemplate, ...],
    source_inputs: dict[str, JsonValue],
    candidate_library_text: str,
) -> ScreenPlan:
    """把已验证候选和模板写入新的不可变筛查目录。"""
    run_dir = config.paths.outputs_dir / "design" / "screens" / run_id
    if run_dir.exists():
        raise DesignError(f"design screen run directory already exists: {run_dir}")
    run_dir.mkdir(parents=True)
    copied_templates = _copy_templates(run_dir=run_dir, templates=templates)
    protocol_path = _write_protocol(run_dir)
    tasks = _write_resfiles_and_tasks(
        config=config,
        run_dir=run_dir,
        candidates=candidates,
        templates=copied_templates,
    )
    snapshot = run_dir / "config.snapshot.txt"
    library_path = run_dir / "candidates.tsv"
    atomic_write_text(snapshot, config.source_snapshot_text)
    atomic_write_text(library_path, candidate_library_text)
    report = config.design.qualification_report
    if report is None:
        raise DesignError("qualified Stage 4 method lacks a qualification report")
    try:
        report_path = report.resolve().relative_to(config.paths.outputs_dir.resolve())
    except ValueError as exc:
        raise DesignError(
            "Stage 4 qualification report must be inside outputs"
        ) from exc
    tool = verify_rosetta_scripts_tool(config.validation.rosetta)
    plan_path = run_dir / SCREEN_PLAN_NAME
    inputs: dict[str, JsonValue] = {
        "config_snapshot": file_record(snapshot, root=run_dir),
        "qualification_report": {
            "path": report_path.as_posix(),
            "sha256": sha256_file(report),
        },
        "candidate_library": file_record(library_path, root=run_dir),
        "protocol": file_record(protocol_path, root=run_dir),
        **source_inputs,
    }
    atomic_write_json(
        plan_path,
        {
            "schema": "vela.design-screen-plan/2",
            "stage": "design_interface_screen",
            "status": "planned",
            "run_id": run_id,
            "planned_at": utc_now(),
            "design_round": design_round,
            "method_id": config.design.method_id,
            "chemistry_id": config.chemistry.chemistry_id,
            "objective": config.design.objective,
            "software": {
                **vela_software_identity(),
                "rosetta_version": tool.version,
                "rosetta_scripts_sha256": tool.executable_sha256,
            },
            "inputs": inputs,
            "parameters": design_parameters(config),
            "candidates": [candidate_record(item) for item in candidates],
            "templates": [
                template_record(item, run_dir=run_dir) for item in copied_templates
            ],
            "tasks": [
                {
                    "task_id": task.task_id,
                    "pair_id": task.pair_id,
                    "state": task.state,
                    "candidate_id": task.candidate.candidate_id,
                    "template_id": task.template.template_id,
                    "seed": task.seed,
                    "resfile": file_record(task.resfile_path, root=run_dir),
                    "wt_context_id": task.wt_context_id,
                    "status": "planned",
                }
                for task in tasks
            ],
        },
    )
    return ScreenPlan(
        run_dir,
        design_round,
        candidates,
        copied_templates,
        tasks,
        protocol_path,
    )


def write_single_screen_plan(
    *,
    config: AppConfig,
    refinement_run_dir: Path,
    run_id: str,
    target_cluster_ids: tuple[str, ...],
) -> ScreenPlan:
    """冻结完整单点库、显式多状态模板和成对 WT 任务矩阵。"""
    validate_run_id(run_id)
    readiness = assess_design_readiness(config)
    if not readiness.screen_ready:
        raise DesignError(
            "Stage 4 is not ready: "
            + "; ".join(issue.code for issue in readiness.issues)
        )
    source = refinement_run_dir.expanduser().resolve()
    source_root = (config.paths.outputs_dir / "validation" / "refinements").resolve()
    try:
        source.relative_to(source_root)
    except ValueError as exc:
        raise DesignError("Stage 4 source is outside Stage 3 refinement runs") from exc
    templates = selected_templates(
        config=config,
        refinement_run_dir=source,
        target_cluster_ids=target_cluster_ids,
    )
    candidates = systematic_single_library(
        chemistry=config.chemistry, settings=config.design
    )
    try:
        source_path = source.relative_to(config.paths.outputs_dir.resolve())
    except ValueError as exc:
        raise DesignError("Stage 4 source record must be inside outputs") from exc
    source_inputs: dict[str, JsonValue] = {
        "stage3_refinement_run": {
            "path": source_path.as_posix(),
            "refinement_manifest_sha256": sha256_file(
                source / "refinement_manifest.json"
            ),
            "analysis_manifest_sha256": sha256_file(
                source / "refinement_analysis" / "analysis_manifest.json"
            ),
            "candidate_review_manifest_sha256": sha256_file(
                source / "candidate_review" / "review_manifest.json"
            ),
        }
    }
    return write_screen_plan(
        config=config,
        run_id=run_id,
        design_round="single",
        candidates=candidates,
        templates=templates,
        source_inputs=source_inputs,
        candidate_library_text=candidate_table(candidates),
    )


def _analysis_candidates(path: Path) -> tuple[dict[str, str], ...]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "candidate_id",
            "mutation_positions",
            "positive_median_paired_dG_separated_delta_REU",
            "positive_worst_paired_dG_separated_delta_REU",
            "candidate_status",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise DesignError("single-screen candidate summary columns are invalid")
        return tuple(dict(row) for row in reader)


def write_combination_screen_plan(
    *, config: AppConfig, single_run_dir: Path, run_id: str
) -> ScreenPlan:
    """从已通过单点门槛的替换建立有限组合并冻结重新评价计划。"""
    validate_run_id(run_id)
    readiness = assess_design_readiness(config)
    if not readiness.screen_ready:
        raise DesignError(
            "Stage 4 is not ready: "
            + "; ".join(issue.code for issue in readiness.issues)
        )
    source = single_run_dir.expanduser().resolve()
    source_root = (config.paths.outputs_dir / "design" / "screens").resolve()
    try:
        source.relative_to(source_root)
        source_path = source.relative_to(config.paths.outputs_dir.resolve())
    except ValueError as exc:
        raise DesignError(
            "single-screen source is outside Stage 4 screen runs"
        ) from exc
    parent = read_screen_plan(config=config, run_dir=source)
    if parent.design_round != "single":
        raise DesignError("combination proposals require a single-mutation screen")
    manifest_path = source / "screen_analysis" / "analysis_manifest.json"
    manifest = read_document(manifest_path, name="single-screen analysis manifest")
    if (
        manifest.get("schema") != "vela.design-screen-analysis-manifest/2"
        or manifest.get("status") != "completed"
        or manifest.get("design_round") != "single"
        or manifest.get("objective") != config.design.objective
    ):
        raise DesignError("single-screen analysis identity is invalid")
    try:
        summary_record = object_mapping(
            manifest.get("candidate_summary"), name="single candidate summary"
        )
    except TypeError as exc:
        raise DesignError("single-screen candidate summary record is invalid") from exc
    summary_path, _ = validate_record(
        root=manifest_path.parent,
        raw=summary_record,
        name="single-screen candidate summary",
    )
    by_id = {item.candidate_id: item for item in parent.candidates}
    ranked: dict[int, list[tuple[float, float, str, str]]] = defaultdict(list)
    for row in _analysis_candidates(summary_path):
        if row["candidate_status"] != "screen_supported":
            continue
        candidate = by_id.get(row["candidate_id"])
        if candidate is None or candidate.mutation_count != 1:
            raise DesignError("single-screen summary references an invalid candidate")
        try:
            worst = float(row["positive_worst_paired_dG_separated_delta_REU"])
            median = float(row["positive_median_paired_dG_separated_delta_REU"])
        except ValueError as exc:
            raise DesignError(
                "single-screen summary contains an invalid score"
            ) from exc
        position = candidate.mutation_positions[0]
        residue = candidate.sequence[position - 1]
        ranked[position].append((worst, median, candidate.candidate_id, residue))
    substitutions = {
        position: tuple(
            item[3]
            for item in sorted(values)[
                : config.design.combination.max_options_per_position
            ]
        )
        for position, values in ranked.items()
    }
    candidates = combination_library(
        chemistry=config.chemistry,
        settings=config.design,
        substitutions=substitutions,
    )
    if not candidates:
        raise DesignError("single-screen evidence produced no legal combinations")
    source_inputs: dict[str, JsonValue] = {
        "parent_single_screen": {
            "path": source_path.as_posix(),
            "screen_plan_sha256": sha256_file(source / SCREEN_PLAN_NAME),
            "analysis_manifest_sha256": sha256_file(manifest_path),
        }
    }
    return write_screen_plan(
        config=config,
        run_id=run_id,
        design_round="combination",
        candidates=candidates,
        templates=parent.templates,
        source_inputs=source_inputs,
        candidate_library_text=candidate_table(candidates),
    )
