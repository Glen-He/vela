"""Rosetta FlexPepDock 工具身份和可复现命令合同。"""

from __future__ import annotations

import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

from vela.core.provenance import sha256_file
from vela.validation.bound_states.assets import PEPTIDE_CHAIN, RECEPTOR_CHAIN
from vela.validation.models import RosettaSettings, ValidationError


@dataclass(frozen=True, slots=True)
class RosettaToolInfo:
    """通过本机检查的 Rosetta 工具身份。"""

    version: str
    executable_sha256: str


@dataclass(frozen=True, slots=True)
class RosettaScriptsToolInfo:
    """通过本机检查的 RosettaScripts 工具身份。"""

    version: str
    executable_sha256: str


def _run_identity_command(command: list[str], *, name: str) -> str:
    result = subprocess.run(command, check=False, capture_output=True, text=True)
    output = (result.stdout + result.stderr).strip()
    if result.returncode != 0 or not output:
        raise ValidationError(f"{name} failed: {output or 'no output'}")
    return output


def run_rosetta_command(
    *,
    command: tuple[str, ...],
    log_path: Path,
    thread_count: int | None = None,
) -> None:
    """执行已构造的 Rosetta 命令并把完整输出写入指定日志。"""
    environment: dict[str, str] | None = None
    if thread_count is not None:
        environment = os.environ.copy()
        for name in (
            "OMP_NUM_THREADS",
            "MKL_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS",
            "NUMEXPR_NUM_THREADS",
        ):
            environment[name] = str(thread_count)
    with log_path.open("w", encoding="utf-8") as log:
        result = subprocess.run(
            command,
            check=False,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
            env=environment,
        )
    if result.returncode != 0:
        raise ValidationError(
            f"Rosetta failed with exit code {result.returncode}: {log_path}"
        )


def single_rosetta_pdb_output(output_dir: Path) -> Path:
    """要求一次单结构 Rosetta 任务恰好生成一个 PDB。"""
    outputs = tuple(output_dir.glob("*.pdb"))
    if len(outputs) != 1:
        raise ValidationError(f"Rosetta must produce exactly one PDB: {output_dir}")
    return outputs[0]


def verify_flexpepdock_tool(settings: RosettaSettings) -> RosettaToolInfo:
    """确认实际二进制、数据库、MPI 和阶段三依赖的命令行能力。"""
    if not settings.executable.is_file() or not os.access(settings.executable, os.X_OK):
        raise ValidationError(
            f"FlexPepDock executable is not runnable: {settings.executable}"
        )
    if not settings.database.is_dir():
        raise ValidationError(f"Rosetta database is missing: {settings.database}")
    if not settings.version_file.is_file():
        raise ValidationError(
            f"Rosetta version declaration is missing: {settings.version_file}"
        )
    source_root = settings.version_file.resolve().parents[2]
    try:
        settings.executable.resolve().relative_to(source_root)
    except ValueError as exc:
        raise ValidationError(
            "FlexPepDock executable and version declaration are from different trees"
        ) from exc
    version_text = settings.version_file.read_text(encoding="utf-8")
    match = re.search(
        r'std::string version\(\)\s*\{\s*return\s+"([^"]+)"', version_text
    )
    if match is None:
        raise ValidationError("FlexPepDock version output is not recognized")
    version = match.group(1)
    if version != settings.expected_version:
        raise ValidationError(
            f"Rosetta version mismatch: expected {settings.expected_version}, got {version}"
        )
    help_text = _run_identity_command(
        [
            str(settings.executable),
            "-database",
            str(settings.database),
            "-help",
        ],
        name="FlexPepDock help check",
    )
    required_options = (
        "fix_disulf",
        "flexpep_prepack",
        "pep_refine",
        "peptide_chain",
        "receptor_chain",
    )
    missing = [option for option in required_options if option not in help_text]
    if missing:
        raise ValidationError(
            "FlexPepDock is missing required options: " + ", ".join(missing)
        )
    return RosettaToolInfo(
        version=version,
        executable_sha256=sha256_file(settings.executable.resolve()),
    )


def verify_rosetta_scripts_tool(
    settings: RosettaSettings,
) -> RosettaScriptsToolInfo:
    """确认 RosettaScripts 与 FlexPepDock 使用同一源码树和数据库。"""
    if not settings.scripts_executable.is_file() or not os.access(
        settings.scripts_executable, os.X_OK
    ):
        raise ValidationError(
            f"RosettaScripts executable is not runnable: {settings.scripts_executable}"
        )
    if not settings.database.is_dir():
        raise ValidationError(f"Rosetta database is missing: {settings.database}")
    if not settings.version_file.is_file():
        raise ValidationError(
            f"Rosetta version declaration is missing: {settings.version_file}"
        )
    source_root = settings.version_file.resolve().parents[2]
    try:
        settings.scripts_executable.resolve().relative_to(source_root)
    except ValueError as exc:
        raise ValidationError(
            "RosettaScripts and version declaration are from different trees"
        ) from exc
    version_text = settings.version_file.read_text(encoding="utf-8")
    match = re.search(
        r'std::string version\(\)\s*\{\s*return\s+"([^"]+)"', version_text
    )
    if match is None or match.group(1) != settings.expected_version:
        raise ValidationError(
            "RosettaScripts version declaration does not match config"
        )
    for mover, required_fields in (
        ("ModifyVariantType", ("ModifyVariantType", "add_type", "residue_selector")),
        (
            "ForceDisulfides",
            ("ForceDisulfides", "disulfides", "scorefxn", "repack"),
        ),
    ):
        mover_schema = _run_identity_command(
            [
                str(settings.scripts_executable),
                "-database",
                str(settings.database),
                "-parser:info",
                mover,
            ],
            name=f"RosettaScripts {mover} schema check",
        )
        if any(field not in mover_schema for field in required_fields):
            raise ValidationError(
                f"RosettaScripts lacks the required {mover} capability"
            )
    return RosettaScriptsToolInfo(
        version=match.group(1),
        executable_sha256=sha256_file(settings.scripts_executable.resolve()),
    )


def build_chemistry_command(
    *,
    settings: RosettaSettings,
    input_path: Path,
    protocol_path: Path,
    disulfide_path: Path,
    output_dir: Path,
    seed: int,
) -> tuple[str, ...]:
    """建立只恢复显式化学变体、不执行局部采样的 RosettaScripts 命令。"""
    if seed < 0:
        raise ValidationError("RosettaScripts seed must not be negative")
    return (
        str(settings.scripts_executable),
        "-database",
        str(settings.database),
        "-s",
        str(input_path),
        "-parser:protocol",
        str(protocol_path),
        "-in:fix_disulf",
        str(disulfide_path),
        "-score:weights",
        settings.score_function,
        "-constant_seed",
        "-jran",
        str(seed),
        "-nstruct",
        "1",
        "-out:path:pdb",
        str(output_dir),
        "-out:path:score",
        str(output_dir),
        "-out:file:scorefile",
        "chemistry.sc",
        "-overwrite",
    )


def _common_command(
    *, settings: RosettaSettings, input_path: Path, disulfide_path: Path
) -> list[str]:
    return [
        str(settings.executable),
        "-database",
        str(settings.database),
        "-s",
        str(input_path),
        "-in:fix_disulf",
        str(disulfide_path),
        "-flexPepDocking:receptor_chain",
        RECEPTOR_CHAIN,
        "-flexPepDocking:peptide_chain",
        PEPTIDE_CHAIN,
        "-score:weights",
        settings.score_function,
        "-ex1",
        "-ex2aro",
        "-use_input_sc",
        "-packing:no_optH",
        "true",
        "-overwrite",
    ]


def build_prepack_command(
    *,
    settings: RosettaSettings,
    input_path: Path,
    disulfide_path: Path,
    output_dir: Path,
    seed: int,
    fixed_histidine_pose_indices: tuple[int, ...] = (),
) -> tuple[str, ...]:
    """建立只重排侧链、不移动复合物骨架的 prepack 命令。"""
    if seed < 0:
        raise ValidationError("FlexPepDock seed must not be negative")
    command = _common_command(
        settings=settings,
        input_path=input_path,
        disulfide_path=disulfide_path,
    )
    if fixed_histidine_pose_indices:
        command.extend(
            [
                "-packing:fix_his_tautomer",
                *(str(index) for index in fixed_histidine_pose_indices),
            ]
        )
    command.extend(
        [
            "-flexPepDocking:flexpep_prepack",
            "-flexPepDocking:flexpep_score_only",
            "-constant_seed",
            "-jran",
            str(seed),
            "-nstruct",
            "1",
            "-out:path:pdb",
            str(output_dir),
            "-out:path:score",
            str(output_dir),
            "-out:file:scorefile",
            "prepack.sc",
        ]
    )
    return tuple(command)


def _refine_command(
    *,
    command: list[str],
    settings: RosettaSettings,
    output_dir: Path,
    seed: int,
    native_path: Path | None = None,
    random_translation_A: float = 0.0,
    random_rotation_degrees: float = 0.0,
    fixed_histidine_pose_indices: tuple[int, ...] = (),
    nstruct: int,
) -> tuple[str, ...]:
    """补齐 MPI 或任务级并行共用的高分辨率精修参数。"""
    if seed < 0:
        raise ValidationError("FlexPepDock seed must not be negative")
    if random_translation_A < 0 or random_rotation_degrees < 0:
        raise ValidationError("FlexPepDock initial perturbations must not be negative")
    if nstruct < 1:
        raise ValidationError("FlexPepDock nstruct must be positive")
    if fixed_histidine_pose_indices:
        command.extend(
            [
                "-packing:fix_his_tautomer",
                *(str(index) for index in fixed_histidine_pose_indices),
            ]
        )
    if native_path is not None:
        command.extend(["-native", str(native_path)])
    if random_translation_A > 0:
        command.extend(
            ["-flexPepDocking:random_trans_start", str(random_translation_A)]
        )
    if random_rotation_degrees > 0:
        command.extend(
            ["-flexPepDocking:random_rot_start", str(random_rotation_degrees)]
        )
    command.extend(["-flexPepDocking:pep_refine"])
    if settings.lowres_preoptimize:
        command.append("-flexPepDocking:lowres_preoptimize")
    command.extend(
        [
            "-flexPepDocking:flexpep_score_only",
            "-constant_seed",
            "-jran",
            str(seed),
            "-nstruct",
            str(nstruct),
            "-out:path:pdb",
            str(output_dir),
            "-out:path:score",
            str(output_dir),
            "-out:file:scorefile",
            "refine.sc",
        ]
    )
    return tuple(command)


def build_refine_command(
    *,
    settings: RosettaSettings,
    input_path: Path,
    disulfide_path: Path,
    output_dir: Path,
    seed: int,
    native_path: Path | None = None,
    random_translation_A: float = 0.0,
    random_rotation_degrees: float = 0.0,
    fixed_histidine_pose_indices: tuple[int, ...] = (),
) -> tuple[str, ...]:
    """建立由外层任务池并行的单进程高分辨率局部精修命令。"""
    return _refine_command(
        command=_common_command(
            settings=settings,
            input_path=input_path,
            disulfide_path=disulfide_path,
        ),
        settings=settings,
        output_dir=output_dir,
        seed=seed,
        native_path=native_path,
        random_translation_A=random_translation_A,
        random_rotation_degrees=random_rotation_degrees,
        fixed_histidine_pose_indices=fixed_histidine_pose_indices,
        nstruct=settings.decoys_per_seed,
    )


def build_topology_refine_command(
    *,
    settings: RosettaSettings,
    input_path: Path,
    disulfide_path: Path,
    output_dir: Path,
    seed: int,
    fixed_histidine_pose_indices: tuple[int, ...] = (),
) -> tuple[str, ...]:
    """建立单结构、无 native 参考的局部拓扑恢复精修命令。"""
    return _refine_command(
        command=_common_command(
            settings=settings,
            input_path=input_path,
            disulfide_path=disulfide_path,
        ),
        settings=settings,
        output_dir=output_dir,
        seed=seed,
        fixed_histidine_pose_indices=fixed_histidine_pose_indices,
        nstruct=1,
    )
