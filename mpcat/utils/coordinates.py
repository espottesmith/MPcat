from enum import Enum
from typing import List, Optional

class ConstraintType(Enum):
    """
    Constraint types corresponding to those defined by schrodinger.infra.mm
    """

    DISTANCE = 1
    ANGLE = 2
    TORSION = 3


class Constraint:
    """
    An object representing a coordinate constraint for a geometry optimization.

    Args:
        constraint_type (ConstraintType): corresponds either to a fixed distance,
            a fixed angle, or a fixed torsion
        atoms (List[int]): List of indices corresponding to atom indices (in
            Pymatgen 0-based indexing). The length of atoms should be 2 for a
            ConstraintType.DISTANCE, 3 for a ConstraintType.ANGLE, and 4
            for a ConstraintType.TORSION
        value (Optional[float]): The value for the constraint. Default is None.
    """

    def __init__(self,
                 constraint_type: ConstraintType,
                 atoms: List[int],
                 value: Optional[float] = None):
        self.constraint_type = constraint_type
        self.atoms = atoms
        self.value = value