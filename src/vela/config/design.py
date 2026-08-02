"""阶段四多状态序列优化配置的严格组装。"""

from collections.abc import Mapping
from pathlib import Path

from vela.config.models import ConfigError
from vela.config.values import (
    assert_keys,
    boolean,
    integer,
    number,
    resolved_path,
    string,
    table,
)
from vela.core.typed_data import object_list
from vela.design.models import (
    CombinationSettings,
    DesignAnalysisSettings,
    DesignSettings,
    FinalistSettings,
    InterfaceScreenSettings,
    IterationSettings,
    SequencePolicy,
)
from vela.discovery.models import UNRESOLVED


def _optional_text(source: Mapping[str, object], key: str, *, path: str) -> str | None:
    value = string(source, key, path=path)
    return None if value == UNRESOLVED else value


def _optional_number(
    source: Mapping[str, object], key: str, *, path: str
) -> float | None:
    if source.get(key) == UNRESOLVED:
        return None
    return number(source, key, path=path)


def _integers(value: object, *, name: str) -> tuple[int, ...]:
    try:
        values = object_list(value, name=name)
    except TypeError as exc:
        raise ConfigError(f"{name} must be an array of integers") from exc
    result: list[int] = []
    for item in values:
        if not isinstance(item, int) or isinstance(item, bool):
            raise ConfigError(f"{name} must be an array of integers")
        result.append(item)
    return tuple(result)


def parse_design(source: Mapping[str, object], *, config_dir: Path) -> DesignSettings:
    """组装 `[design]` 及其序列、筛查、组合和分析子节。"""
    section = table(source, "design", path="")
    required = {
        "method_id",
        "qualification_status",
        "qualification_report",
        "qualification_report_sha256",
        "objective",
        "seeds",
        "sequence",
        "screen",
        "combination",
        "iteration",
        "analysis",
        "finalists",
    }
    assert_keys(section, allowed=required, required=required, path="design")

    sequence = table(section, "sequence", path="design")
    sequence_keys = {
        "mutable_positions",
        "allowed_amino_acids",
        "candidate_histidine_state",
    }
    assert_keys(
        sequence,
        allowed=sequence_keys,
        required=sequence_keys,
        path="design.sequence",
    )

    screen = table(section, "screen", path="design")
    screen_keys = {
        "parallel_tasks",
        "neighbor_distance_A",
        "score_function",
        "ranking_score",
        "pack_separated",
    }
    assert_keys(screen, allowed=screen_keys, required=screen_keys, path="design.screen")

    combination = table(section, "combination", path="design")
    combination_keys = {
        "min_mutations",
        "max_mutations",
        "max_options_per_position",
        "max_candidates",
    }
    assert_keys(
        combination,
        allowed=combination_keys,
        required=combination_keys,
        path="design.combination",
    )

    iteration = table(section, "iteration", path="design")
    iteration_keys = {
        "max_parents",
        "max_total_mutations",
        "max_candidates",
    }
    assert_keys(
        iteration,
        allowed=iteration_keys,
        required=iteration_keys,
        path="design.iteration",
    )

    analysis = table(section, "analysis", path="design")
    analysis_keys = {
        "calibrated",
        "max_positive_median_delta",
        "max_positive_worst_delta",
    }
    assert_keys(
        analysis,
        allowed=analysis_keys,
        required=analysis_keys,
        path="design.analysis",
    )

    finalists = table(section, "finalists", path="design")
    finalist_keys = {
        "parallel_tasks",
        "max_candidates",
        "max_md_candidates",
        "seeds",
        "ranking_score",
        "interface_score",
        "calibrated",
        "min_passed_decoy_fraction",
        "min_successful_seeds",
        "max_positive_median_ranking_delta",
        "max_positive_worst_ranking_delta",
        "max_positive_median_interface_delta",
        "max_positive_worst_interface_delta",
    }
    assert_keys(
        finalists,
        allowed=finalist_keys,
        required=finalist_keys,
        path="design.finalists",
    )

    report = _optional_text(section, "qualification_report", path="design")
    return DesignSettings(
        method_id=_optional_text(section, "method_id", path="design"),
        qualification_status=string(section, "qualification_status", path="design"),
        qualification_report=(
            None if report is None else resolved_path(report, config_dir=config_dir)
        ),
        qualification_report_sha256=_optional_text(
            section, "qualification_report_sha256", path="design"
        ),
        objective=_optional_text(section, "objective", path="design"),
        seeds=_integers(section["seeds"], name="design.seeds"),
        sequence=SequencePolicy(
            mutable_positions=_integers(
                sequence["mutable_positions"],
                name="design.sequence.mutable_positions",
            ),
            allowed_amino_acids=string(
                sequence, "allowed_amino_acids", path="design.sequence"
            ),
            candidate_histidine_state=string(
                sequence, "candidate_histidine_state", path="design.sequence"
            ),
        ),
        screen=InterfaceScreenSettings(
            parallel_tasks=integer(screen, "parallel_tasks", path="design.screen"),
            neighbor_distance_A=number(
                screen, "neighbor_distance_A", path="design.screen"
            ),
            score_function=string(screen, "score_function", path="design.screen"),
            ranking_score=string(screen, "ranking_score", path="design.screen"),
            pack_separated=boolean(screen, "pack_separated", path="design.screen"),
        ),
        combination=CombinationSettings(
            min_mutations=integer(
                combination, "min_mutations", path="design.combination"
            ),
            max_mutations=integer(
                combination, "max_mutations", path="design.combination"
            ),
            max_options_per_position=integer(
                combination,
                "max_options_per_position",
                path="design.combination",
            ),
            max_candidates=integer(
                combination, "max_candidates", path="design.combination"
            ),
        ),
        iteration=IterationSettings(
            max_parents=integer(iteration, "max_parents", path="design.iteration"),
            max_total_mutations=integer(
                iteration,
                "max_total_mutations",
                path="design.iteration",
            ),
            max_candidates=integer(
                iteration, "max_candidates", path="design.iteration"
            ),
        ),
        analysis=DesignAnalysisSettings(
            calibrated=boolean(analysis, "calibrated", path="design.analysis"),
            max_positive_median_delta=_optional_number(
                analysis, "max_positive_median_delta", path="design.analysis"
            ),
            max_positive_worst_delta=_optional_number(
                analysis, "max_positive_worst_delta", path="design.analysis"
            ),
        ),
        finalists=FinalistSettings(
            parallel_tasks=integer(
                finalists, "parallel_tasks", path="design.finalists"
            ),
            max_candidates=integer(
                finalists, "max_candidates", path="design.finalists"
            ),
            max_md_candidates=integer(
                finalists, "max_md_candidates", path="design.finalists"
            ),
            seeds=_integers(finalists["seeds"], name="design.finalists.seeds"),
            ranking_score=string(finalists, "ranking_score", path="design.finalists"),
            interface_score=string(
                finalists, "interface_score", path="design.finalists"
            ),
            calibrated=boolean(finalists, "calibrated", path="design.finalists"),
            min_passed_decoy_fraction=_optional_number(
                finalists,
                "min_passed_decoy_fraction",
                path="design.finalists",
            ),
            min_successful_seeds=(
                None
                if finalists.get("min_successful_seeds") == UNRESOLVED
                else integer(
                    finalists,
                    "min_successful_seeds",
                    path="design.finalists",
                )
            ),
            max_positive_median_ranking_delta=_optional_number(
                finalists,
                "max_positive_median_ranking_delta",
                path="design.finalists",
            ),
            max_positive_worst_ranking_delta=_optional_number(
                finalists,
                "max_positive_worst_ranking_delta",
                path="design.finalists",
            ),
            max_positive_median_interface_delta=_optional_number(
                finalists,
                "max_positive_median_interface_delta",
                path="design.finalists",
            ),
            max_positive_worst_interface_delta=_optional_number(
                finalists,
                "max_positive_worst_interface_delta",
                path="design.finalists",
            ),
        ),
    )
