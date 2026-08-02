"""阶段二标准环肽阳性控制的唯一结构和化学定义。"""

from vela.config import AppConfig
from vela.discovery.models import DiscoveryError
from vela.preparation.chemistry import ChemistryDefinition, HistidineState
from vela.validation.models import BoundStateDefinition


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
    """建立只用于粗粒化采样和拓扑校准的控制肽化学身份。"""
    sequence = state.ligand_sequence
    if sequence is None:
        raise DiscoveryError("qualification control ligand sequence is unresolved")
    return ChemistryDefinition(
        ligand_id=state.ligand_id.lower(),
        chemistry_id=f"qualification-{state.state_id}",
        sequence=sequence,
        chirality="L",
        disulfide_bonds=state.disulfide_bonds,
        n_terminus="not_assessed_by_topology_calibration",
        c_terminus="not_assessed_by_topology_calibration",
        target_ph=None,
        net_charge=None,
        histidines=tuple(
            HistidineState(item.position, item.state) for item in state.histidines
        ),
        other_modifications_status="none",
        other_modifications=(),
        decision_sources=(f"experimental control {state.state_id}",),
    )
