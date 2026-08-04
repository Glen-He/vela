"""阶段二 CABS-dock 外部程序合同、命令和原生二硫键校验。"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

from vela.core.provenance import JsonValue, sha256_file
from vela.discovery.models import CabsDockSettings, DiscoveryError, DiscoveryTask
from vela.preparation.chemistry import ChemistryDefinition

CRITICAL_CABS_SOURCES = (
    "CABS/analysis/restraints.py",
    "CABS/core/cabs.py",
    "CABS/core/job.py",
    "CABS/core/trajectory.py",
    "CABS/data/data0.dat",
    "CABS/io/config.json",
    "CABS/structures/atom.py",
    "CABS/utils/filter.py",
    "CABS/utils/utils.py",
)
CABS_TASK_RESULT_SCHEMA = "vela.cabsdock-task-result/5"


def cabsdock_source_records(
    settings: CabsDockSettings,
) -> dict[str, dict[str, JsonValue]]:
    """冻结会决定 TRAF 解码和 Top-1000 能量映射的上游源码。"""
    records: dict[str, dict[str, JsonValue]] = {}
    for relative in CRITICAL_CABS_SOURCES:
        path = settings.source_dir / relative
        if not path.is_file():
            raise DiscoveryError(f"critical CABS source is missing: {path}")
        records[relative] = {
            "path": relative,
            "sha256": sha256_file(path),
        }
    return records


def _tool_version(executable: Path) -> str:
    result = subprocess.run(
        [str(executable), "--version"],
        check=False,
        capture_output=True,
        text=True,
    )
    output = (result.stdout + result.stderr).strip()
    if result.returncode != 0 or not output:
        raise DiscoveryError(f"CABS-dock version check failed: {output}")
    return output


def verify_cabsdock_tool(settings: CabsDockSettings) -> str:
    executable = settings.executable
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise DiscoveryError(f"CABS-dock executable is not runnable: {executable}")
    version_text = _tool_version(executable)
    help_result = subprocess.run(
        [str(executable), "--help"],
        check=False,
        capture_output=True,
        text=True,
    )
    help_text = help_result.stdout + help_result.stderr
    if help_result.returncode != 0 or "--ca-rest-add" not in help_text:
        raise DiscoveryError(
            "CABS-dock does not expose the required --ca-rest-add option"
        )
    if not settings.source_dir.is_dir():
        raise DiscoveryError(f"CABS source directory is missing: {settings.source_dir}")
    cabsdock_source_records(settings)
    python_executable = executable.parent / "python"
    if not python_executable.is_file():
        raise DiscoveryError(f"CABS Python executable is missing: {python_executable}")
    import_result = subprocess.run(
        [
            str(python_executable),
            "-c",
            (
                "from pathlib import Path; import CABS; "
                "print(Path(CABS.__file__).resolve().parent.parent)"
            ),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    imported_source = import_result.stdout.strip()
    if (
        import_result.returncode != 0
        or not imported_source
        or Path(imported_source).resolve() != settings.source_dir.resolve()
    ):
        raise DiscoveryError(
            "CABS executable does not import the declared source directory: "
            f"expected {settings.source_dir}, got {imported_source or 'unavailable'}"
        )
    revision_result = subprocess.run(
        ["git", "-C", str(settings.source_dir), "rev-parse", "HEAD"],
        check=False,
        capture_output=True,
        text=True,
    )
    revision = revision_result.stdout.strip()
    if revision_result.returncode != 0 or revision != settings.source_revision:
        raise DiscoveryError(
            "CABS source revision mismatch: "
            f"expected {settings.source_revision}, got {revision or 'unavailable'}"
        )
    diff_result = subprocess.run(
        [
            "git",
            "-C",
            str(settings.source_dir),
            "diff",
            "--quiet",
            "HEAD",
            "--",
            *CRITICAL_CABS_SOURCES,
        ],
        check=False,
    )
    if diff_result.returncode == 1:
        raise DiscoveryError(
            "critical CABS source files contain uncommitted changes; restore or commit "
            "the declared source revision before running Stage 2"
        )
    if diff_result.returncode != 0:
        raise DiscoveryError("CABS critical-source worktree check failed")
    return version_text


def build_cabsdock_command(
    *,
    task: DiscoveryTask,
    settings: CabsDockSettings,
    chemistry: ChemistryDefinition,
    secondary_structure: str,
    task_dir: Path,
) -> tuple[str, ...]:
    """建立无已知位点约束、关闭全原子重建的 CABS-dock 命令。"""
    if len(secondary_structure) != len(chemistry.sequence):
        raise DiscoveryError(
            "CABS-dock peptide secondary structure length does not match the ligand"
        )
    command = [
        str(settings.executable),
        "-i",
        str(task.receptor_path),
        "-p",
        f"{chemistry.sequence}:{secondary_structure}",
    ]
    for bond in chemistry.disulfide_bonds:
        command.extend(
            [
                "--ca-rest-add",
                f"{bond.first}:PEP1",
                f"{bond.second}:PEP1",
                str(settings.disulfide_ca_restraint_distance_A),
                str(settings.disulfide_ca_restraint_weight),
            ]
        )
    command.extend(
        [
            "-g",
            "rigid",
            str(settings.protein_restraint_gap),
            str(settings.protein_restraint_min_A),
            str(settings.protein_restraint_max_A),
            "-a",
            str(settings.mc_annealing),
            "-y",
            str(settings.mc_cycles),
            "-s",
            str(settings.mc_steps),
            "-r",
            str(settings.replicas),
            "-D",
            str(settings.replicas_dtemp),
            "-t",
            str(settings.temperature_initial),
            str(settings.temperature_final),
            "-b",
            str(settings.binding_interactions),
            "-n",
            str(settings.filtering_count),
            "--filtering-mode",
            "each",
            "-k",
            str(settings.clustering_medoids),
            "--clustering-iterations",
            str(settings.clustering_iterations),
            "-z",
            str(task.seed),
            "-A",
            "N",
            "-o",
            "FMS",
            "-S",
            "-C",
            "--json-output",
            "--restraints-output",
            "--no-progress-bar",
            "-w",
            str(task_dir),
        ]
    )
    return tuple(command)


def verify_disulfide_ca_restraint(
    *,
    task_dir: Path,
    chemistry: ChemistryDefinition,
    settings: CabsDockSettings,
) -> None:
    """确认任务实际使用了项目声明的 C-alpha 环拓扑约束。"""
    config_text = (task_dir / "config.ini").read_text(encoding="utf-8")
    restraints_text = (task_dir / "output_data" / "restraints.txt").read_text(
        encoding="utf-8"
    )
    for bond in chemistry.disulfide_bonds:
        endpoints = f"{bond.first}:PEP1 {bond.second}:PEP1"
        configured = (
            f"{endpoints} {settings.disulfide_ca_restraint_distance_A} "
            f"{settings.disulfide_ca_restraint_weight}"
        )
        if configured not in config_text:
            raise DiscoveryError(
                f"CABS-dock did not record the disulfide CA restraint: {endpoints}"
            )
        restraint = (
            f"{endpoints} {settings.disulfide_ca_restraint_distance_A:.4f} "
            f"1.0000 {settings.disulfide_ca_restraint_weight:.2f}"
        )
        if restraint not in restraints_text:
            raise DiscoveryError(
                f"CABS-dock did not generate the disulfide CA restraint: {endpoints}"
            )


def cabsdock_medoid_paths(task_dir: Path) -> tuple[Path, ...]:
    """按 CABS-dock 的零基编号返回独立 medoid 文件。"""
    output_dir = task_dir / "output_pdbs"
    numbered: list[tuple[int, Path]] = []
    for path in output_dir.glob("model_*.pdb"):
        suffix = path.stem.removeprefix("model_")
        if not suffix.isdigit():
            raise DiscoveryError(f"invalid CABS-dock medoid filename: {path}")
        numbered.append((int(suffix), path))
    if numbered:
        numbered.sort(key=lambda item: item[0])
        indices = tuple(index for index, _ in numbered)
        if indices != tuple(range(len(numbered))):
            raise DiscoveryError(
                f"CABS-dock medoid indices are not contiguous: {output_dir}"
            )
        if (output_dir / "model.pdb").exists():
            raise DiscoveryError(
                f"CABS-dock output contains ambiguous medoid files: {output_dir}"
            )
        return tuple(path for _, path in numbered)
    singleton = output_dir / "model.pdb"
    if singleton.is_file():
        return (singleton,)
    raise DiscoveryError(f"CABS-dock medoid output is missing: {output_dir}")


def cabsdock_archive_path(task_dir: Path) -> Path:
    """返回单任务唯一的原始 CABS 归档。"""
    archives = sorted(task_dir.glob("*.cbs"))
    if len(archives) != 1:
        raise DiscoveryError(
            f"CABS-dock task must contain one .cbs archive: {task_dir}"
        )
    return archives[0]


def cabsdock_output_records(task_dir: Path) -> dict[str, dict[str, JsonValue]]:
    paths = {
        "config": task_dir / "config.ini",
        "restraints": task_dir / "output_data" / "restraints.txt",
        "filtered_models": task_dir / "output_pdbs" / "top1000.pdb",
        "log": task_dir / "cabsdock.log",
    }
    for index, path in enumerate(cabsdock_medoid_paths(task_dir), 1):
        paths[f"medoid_{index:02d}"] = path
    paths["cabs_archive"] = cabsdock_archive_path(task_dir)
    paths["trajectory_audit"] = task_dir / "trajectory_audit.json"
    records: dict[str, dict[str, JsonValue]] = {}
    for name, path in paths.items():
        if not path.is_file():
            raise DiscoveryError(f"CABS-dock output is missing: {path}")
        records[name] = {
            "path": path.relative_to(task_dir).as_posix(),
            "sha256": sha256_file(path),
        }
    return records
