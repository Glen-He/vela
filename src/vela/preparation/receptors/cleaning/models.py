"""受体基础制备的结果和替代构象决定。"""

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class PreparationResult:
    """一个受体基础模型的制备统计。"""

    receptor_id: str
    pdb_id: str
    output_path: Path
    source_sha256: str
    output_sha256: str
    residue_count: int
    atom_count: int
    removed_water_residues: int
    removed_nonpolymer_residues: int
    altloc_residues_resolved: int
    altloc_atoms_removed: int


@dataclass(frozen=True, slots=True)
class AltlocDecision:
    """一个残基的替代构象选择记录。"""

    receptor_id: str
    residue_id: str
    residue_name: str
    available_altlocs: str
    selected_altloc: str
    removed_atom_count: int
