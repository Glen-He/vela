from pathlib import Path

import pytest

from vela.core.errors import VelaError
from vela.core.provenance import (
    JsonValue,
    atomic_write_json,
    is_current_vela_software,
    is_vela_software_identity,
    vela_software_identity,
)
from vela.design.models import DesignError
from vela.validation.models import ValidationError
from vela.validation.records import (
    ResumedResult,
    file_record,
    resume_completed_result,
)


def _write_completed_result(directory: Path) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    output = directory / "output.pdb"
    output.write_text("MODEL\nENDMDL\n", encoding="utf-8")
    result_path = directory / "task_result.json"
    atomic_write_json(
        result_path,
        {
            "schema": "vela.test-task-result/1",
            "status": "completed",
            "task_id": "task_001",
            "pair_id": "pair_001",
            "test_plan_sha256": "plan-hash",
            "output": file_record(output, root=directory),
        },
    )
    return result_path


def _resume(
    directory: Path, *, stale_error: type[VelaError] = ValidationError
) -> ResumedResult | None:
    return resume_completed_result(
        directory=directory,
        filename="task_result.json",
        document_name="test task result",
        schema="vela.test-task-result/1",
        identity={"task_id": "task_001", "pair_id": "pair_001"},
        plan_hash_key="test_plan_sha256",
        plan_hash="plan-hash",
        records={"output": "test output"},
        stale_label="test task",
        stale_error=stale_error,
    )


def test_resume_completed_result_returns_none_when_checkpoint_is_absent(
    tmp_path: Path,
) -> None:
    assert _resume(tmp_path) is None


def test_atomic_json_rejects_non_finite_numbers(tmp_path: Path) -> None:
    target = tmp_path / "invalid.json"

    with pytest.raises(ValueError, match="Out of range float values"):
        atomic_write_json(target, {"metric": float("nan")})

    assert not target.exists()


def test_software_identity_detects_a_different_source_tree() -> None:
    identity = vela_software_identity()

    assert is_vela_software_identity(identity)
    assert is_vela_software_identity({**identity, "vela_source_sha256": "0" * 64})
    assert is_current_vela_software(identity)
    assert not is_current_vela_software({**identity, "vela_source_sha256": "0" * 64})


@pytest.mark.parametrize(
    "identity",
    [
        {},
        {"vela_version": "0.1.0", "vela_source_sha256": "short"},
        {"vela_version": "", "vela_source_sha256": "0" * 64},
        {"vela_version": "0.1.0", "vela_source_sha256": "G" * 64},
    ],
)
def test_software_identity_rejects_invalid_historical_records(
    identity: dict[str, str],
) -> None:
    assert not is_vela_software_identity(identity)


def test_resume_completed_result_validates_identity_and_files(tmp_path: Path) -> None:
    result_path = _write_completed_result(tmp_path)

    resumed = _resume(tmp_path)

    assert resumed is not None
    assert resumed.path == result_path
    assert resumed.document["task_id"] == "task_001"
    assert resumed.files == {"output": tmp_path / "output.pdb"}


@pytest.mark.parametrize(
    ("field", "replacement"),
    [
        ("schema", "vela.wrong/1"),
        ("status", "running"),
        ("task_id", "task_002"),
        ("pair_id", "pair_002"),
        ("test_plan_sha256", "other-plan-hash"),
    ],
)
def test_resume_completed_result_rejects_stale_checkpoint(
    tmp_path: Path, field: str, replacement: JsonValue
) -> None:
    result_path = _write_completed_result(tmp_path)
    document: dict[str, JsonValue] = {
        "schema": "vela.test-task-result/1",
        "status": "completed",
        "task_id": "task_001",
        "pair_id": "pair_001",
        "test_plan_sha256": "plan-hash",
        "output": file_record(tmp_path / "output.pdb", root=tmp_path),
    }
    document[field] = replacement
    atomic_write_json(result_path, document)

    with pytest.raises(ValidationError, match="stale test task: task_001"):
        _resume(tmp_path)


def test_resume_completed_result_preserves_stage_error_type(tmp_path: Path) -> None:
    result_path = _write_completed_result(tmp_path)
    atomic_write_json(
        result_path,
        {
            "schema": "vela.test-task-result/1",
            "status": "completed",
            "task_id": "task_001",
            "pair_id": "pair_001",
            "test_plan_sha256": "stale-plan-hash",
            "output": file_record(tmp_path / "output.pdb", root=tmp_path),
        },
    )

    with pytest.raises(DesignError, match="stale test task: task_001"):
        _resume(tmp_path, stale_error=DesignError)


def test_resume_completed_result_rejects_changed_output(tmp_path: Path) -> None:
    _write_completed_result(tmp_path)
    (tmp_path / "output.pdb").write_text("CHANGED\n", encoding="utf-8")

    with pytest.raises(ValidationError, match="test output hash mismatch"):
        _resume(tmp_path)
