"""阶段二标准环肽阳性控制的唯一结构和化学定义。"""

from vela.config import AppConfig
from vela.discovery.models import DiscoveryError
from vela.preparation.chemistry import ChemistryDefinition
from vela.validation.bound_states.assets import standard_peptide_chemistry
from vela.validation.models import BoundStateDefinition, ValidationError


def control_bound_state(config: AppConfig) -> BoundStateDefinition:
    """返回与声明靶标一致的标准环肽实验复合物。"""
    qualification = config.discovery.qualification
    for state in config.validation.bound_states:
        if state.state_id != qualification.control_bound_state_id:
            continue
        if state.local_control_kind != "standard_cyclic_peptide":
            raise DiscoveryError(
                "qualification control must be a standard cyclic peptide"
            )
        receptor = next(
            (
                item
                for item in config.receptors
                if item.receptor_id == state.receptor_id
            ),
            None,
        )
        if receptor is None or receptor.target != qualification.control_target_id:
            raise DiscoveryError(
                "qualification control target differs from its receptor definition"
            )
        return state
    raise DiscoveryError(
        "qualification control bound state is unknown: "
        f"{qualification.control_bound_state_id}"
    )


def control_chemistry(state: BoundStateDefinition) -> ChemistryDefinition:
    """建立贯穿粗粒化采样、拓扑恢复和局部资格的控制肽化学身份。"""
    try:
        return standard_peptide_chemistry(state)
    except ValidationError as exc:
        raise DiscoveryError(str(exc)) from exc
