from CABS.structures.atom import Atoms

class Header:
    model: int
    replica: int

class Coordinates:
    def reshape(self, *shape: int) -> Coordinates: ...
    def __getitem__(self, index: int) -> object: ...

class Trajectory:
    template: Atoms
    coordinates: Coordinates
    headers: list[Header]

    @classmethod
    def read_trajectory(cls, traf: str, seq: str) -> Trajectory: ...
