"""多受体 receptor-only 基础制备编排。"""

from pathlib import Path

from vela.preparation.receptors.cleaning.models import (
    AltlocDecision,
    PreparationResult,
)
from vela.preparation.receptors.cleaning.reports import write_preparation_reports
from vela.preparation.receptors.cleaning.structure import prepare_structure
from vela.preparation.receptors.models import (
    ReceptorDefinition,
    ReceptorPreparationConfig,
)


def prepare_receptors(
    *,
    definitions: tuple[ReceptorDefinition, ...],
    settings: ReceptorPreparationConfig,
    data_dir: Path,
) -> tuple[PreparationResult, ...]:
    """制备所有 prepare=true 的受体并记录父文件和全部处理决定。"""
    raw_dir = data_dir / "receptors" / "raw"
    prepared_dir = data_dir / "receptors" / "prepared"
    results: list[PreparationResult] = []
    decisions: list[AltlocDecision] = []
    for definition in definitions:
        if not definition.prepare:
            continue
        result, entry_decisions = prepare_structure(
            definition=definition,
            altloc_settings=settings.altloc,
            raw_dir=raw_dir,
            prepared_dir=prepared_dir,
        )
        results.append(result)
        decisions.extend(entry_decisions)
    write_preparation_reports(
        results=results,
        decisions=decisions,
        definitions=definitions,
        settings=settings,
        data_dir=data_dir,
    )
    return tuple(results)
