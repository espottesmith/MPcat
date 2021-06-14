# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import logging

from monty.io import zopen
from monty.json import MSONable

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN, metal_edge_extender
from pymatgen.io.qchem.inputs import QCInput

from pymatgen.io.qchem.utils import lower_and_check_unique

# Classes for reading/manipulating/writing input files for use with pyGSM.

__author__ = "Evan Spotte-Smith"
__copyright__ = "Copyright 2021"
__version__ = "0.1"
__email__ = "espottesmith@gmail.com"

logger = logging.getLogger(__name__)


class QCTemplate(QCInput):
    """
    An object representing an incomplete Q-Chem input file used as a template
    for calculations (especially force/gradient calculations) in the Growing
    String Method (GSM), as implemented in pyGSM.

    Args:
        rem (dict): A dictionary of all the input parameters for the rem section
            of Q-Chem input file.
            Ex: rem = {'method': 'rimp2', 'basis': '6-31*G++' ... }
        pcm (dict): A dictionary of values relating to the polarizable continuum
            model (PCM). Note that, if a pcm dict is provided, then a "solvent"
            dict (described below) should also be provided, but a "smx" dict
            (also described below, and used for SMX methods like SMD) should not
            be provided.
            Ex: pcm = {"theory": "ief-pcm", "radii": "uff"}
        solvent (dict): A dictionary of values related to the solvent used in
            the PCM method. Note that this section should be provided if PCM
            is to be used (see above), but should not be provided for use with
            the SMD or other SMX methods.
            Ex: solvent = {"dielectric": 78.39}
        smx (dict): A dictionary of values for use with the SMX methods (like
            SMD). If this section is provided, then neither the pcm nor the
            solvent sections should be provided.
    """

    def __init__(self, rem, pcm=None, solvent=None, smx=None):
        self.rem = lower_and_check_unique(rem)
        self.pcm = lower_and_check_unique(pcm)
        self.solvent = lower_and_check_unique(solvent)
        self.smx = lower_and_check_unique(smx)

        # Only force jobs should be used with GSM
        self.rem["job_type"] = "force"

        if "basis" not in self.rem:
            raise ValueError("The rem dictionary must contain a 'basis' entry")
        if "method" not in self.rem:
            if "exchange" not in self.rem:
                raise ValueError(
                    "The rem dictionary must contain either a 'method' entry or an 'exchange' entry"
                )

    def __str__(self):
        combined_list = list()
        # rem section
        combined_list.append(self.rem_template(self.rem))
        combined_list.append("")
        # pcm section
        if self.pcm:
            combined_list.append(self.pcm_template(self.pcm))
            combined_list.append("")
        # solvent section
        if self.solvent:
            combined_list.append(self.solvent_template(self.solvent))
            combined_list.append("")
        if self.smx:
            combined_list.append(self.smx_template(self.smx))
            combined_list.append("")

        # Finally, add beginning of molecule section
        combined_list.append("$molecule")

        return '\n'.join(combined_list) + "\n"

    @classmethod
    def from_string(cls, string):
        sections = cls.find_sections(string)
        rem = cls.read_rem(string)
        # only rem is necessary everything else is checked
        pcm = None
        solvent = None
        smx = None
        if "pcm" in sections:
            pcm = cls.read_pcm(string)
        if "solvent" in sections:
            solvent = cls.read_solvent(string)
        if "smx" in sections:
            smx = cls.read_smx(string)

        return cls(rem, pcm=pcm, solvent=solvent, smx=smx)

    @staticmethod
    def from_file(filename):
        with zopen(filename, 'rt') as f:
            return QCTemplate.from_string(f.read())


class GSMIsomerInput(MSONable):
    """
    An object representing an isomers input file for single-ended growing-string
    method (GSM) calculations.

    NOTE: All indices here are zero-based, reflecting the indexing scheme used
        in the pymatgen Molecule class and in Python generally.

    Args:
        molecule (pymatgen.Molecule object): Molecule that will be subjected to
            analysis with GSM
        bonds_formed (list of tuples): list of bond formation driving
            coordinates
            Ex: bonds_formed = [(1, 3), (2, 10)] indicates that atoms 1 and 3
            should bond (move closer together), as should atoms 2 and 10.
        bonds_broken (list of tuples): list of bond cleavage driving coordinates
            Ex: bonds_broken = [(1, 5)] indicates that atoms 1 and 5, initially
            bonded, should break apart.
        angles (list of tuples): list of angle driving coordinates
            Ex: angles = [(1, 4, 5)] indicates that the angle between atoms
            1, 4, and 5 should be varied.
        torsions (list of tuples): list of dihedral angle driving coordinates
            Ex: torsions = [(3, 6, 7, 10)] indicates that the dihedral angle
            between atoms 3, 6, 7, and 10 should be varied.
        out_of_planes (list of tuples): list of out-of-plane bend driving
            coordinates.
            Ex: out_of_planes = [(1, 3, 8, 9), (2, 4, 8, 11)] indicates that the
            out-of-plane bends for atoms 1, 3, 8, and 9, and between atoms 2, 4,
            8, and 11, should be varied.
        use_graph (bool): if True (default False), verify that the provided
            coordinates are valid, based on the current connectivity of the
            molecule. Even if this is False, some checks will be made.

    """

    def __init__(self, molecule=None, bonds_formed=None, bonds_broken=None,
                 angles=None, torsions=None, out_of_planes=None, use_graph=False):
        self.molecule = molecule
        self.bonds_formed = bonds_formed or list()
        self.bonds_broken = bonds_broken or list()
        self.angles = angles or list()
        self.torsions = torsions or list()
        self.out_of_planes = out_of_planes or list()
        self.use_graph = use_graph

        # First, check that there are not too many coordinates given
        # We do not allow more than 4 driving coordinates
        num_coords = 0
        all_coords = list()

        if self.bonds_formed is not None:
            num_coords += len(self.bonds_formed)
            all_coords += self.bonds_formed
        if self.bonds_broken is not None:
            num_coords += len(self.bonds_broken)
            all_coords += self.bonds_broken
        if self.angles is not None:
            num_coords += len(self.angles)
            all_coords += self.angles
        if self.torsions is not None:
            num_coords += len(self.torsions)
            all_coords += self.torsions
        if self.out_of_planes is not None:
            num_coords += len(self.out_of_planes)
            all_coords += self.out_of_planes

        if num_coords > 4:
            raise ValueError("Too many driving coordinates given! At most 4 "
                             "driving coordinates may be used with "
                             "GSMIsomerInput.")
        else:
            self.num_coords = num_coords
            self.all_coords = all_coords

        # Verify that all coordinates are tuples including only valid indices
        for coord in all_coords:
                for index in coord:
                    if isinstance(index, int):
                        if self.molecule is not None:
                            if index >= len(self.molecule):
                                raise ValueError("Invalid index given for coordinate {}!".format(coord))
                    else:
                        raise ValueError("Non-integer index given for coordinate {}!".format(coord))

        # Verify that coordinates are of appropriate length for their type
        for coord in self.bonds_broken + self.bonds_formed:
            if len(set(coord)) != 2:
                raise ValueError("All bond coordinates should involve 2 indices!")
        for coord in self.angles:
            if len(set(coord)) != 3:
                raise ValueError("All angle coordinates should involve 3 indices!")
        for coord in self.torsions + self.out_of_planes:
            if len(set(coord)) != 4:
                raise ValueError("All angle coordinates should involve 4 indices!")

        # If allowed, use graph methods to verify the coordinates
        if self.use_graph and self.molecule is not None:
            self.molecule_graph = MoleculeGraph.with_local_env_strategy(self.molecule,
                                                                        OpenBabelNN())
            self.molecule_graph = metal_edge_extender(self.molecule_graph)

            self._verify_two_atom_coords()
            self._verify_three_atom_coords()
            self._verify_four_atom_coords()

        else:
            self.molecule_graph = None

    def _verify_two_atom_coords(self):
        for atom_1, atom_2 in self.bonds_formed:
            if atom_2 in self.molecule_graph.graph[atom_1] or atom_1 in self.molecule_graph.graph[atom_2]:
                raise ValueError("Requested coordinate in bonds_formed already exists!")

        for atom_1, atom_2 in self.bonds_broken:
            if atom_2 not in self.molecule_graph.graph[atom_1] and atom_1 not in self.molecule_graph.graph[atom_2]:
                raise ValueError("Requested coordinate in bonds_broken does not exist!")

    def _verify_three_atom_coords(self):
        for atom_1, atom_2, atom_3 in self.angles:
            first_bond = True
            second_bond = True

            if atom_2 not in self.molecule_graph.graph[atom_1] and atom_1 not in self.molecule_graph.graph[atom_2]:
                first_bond = False
            if atom_2 not in self.molecule_graph.graph[atom_3] and atom_3 not in self.molecule_graph.graph[atom_2]:
                second_bond = False

            if not (first_bond and second_bond):
                raise ValueError("Atoms for requested angle coordinate are not all connected!")

    def _verify_four_atom_coords(self):
        for atom_1, atom_2, atom_3, atom_4 in self.torsions + self.out_of_planes:
            first_bond = True
            second_bond = True
            third_bond = True

            if atom_2 not in self.molecule_graph.graph[atom_1] and atom_1 not in self.molecule_graph.graph[atom_2]:
                first_bond = False
            if atom_2 not in self.molecule_graph.graph[atom_3] and atom_3 not in self.molecule_graph.graph[atom_2]:
                second_bond = False
            if atom_4 not in self.molecule_graph.graph[atom_3] and atom_3 not in self.molecule_graph.graph[atom_4]:
                third_bond = False

            if not (first_bond and second_bond and third_bond):
                raise ValueError("Atoms for requested torsion or out-of-plane coordinate are not all connected!")

    def __str__(self):
        combined = list()

        if self.bonds_formed is not None:
            for atom_1, atom_2 in self.bonds_formed:
                combined.append("ADD {} {}".format(atom_1 + 1, atom_2 + 1))
        if self.bonds_broken is not None:
            for atom_1, atom_2 in self.bonds_broken:
                combined.append("BREAK {} {}".format(atom_1 + 1, atom_2 + 1))
        if self.angles is not None:
            for atom_1, atom_2, atom_3 in self.angles:
                combined.append("ANGLE {} {} {}".format(atom_1 + 1,
                                                        atom_2 + 1,
                                                        atom_3 + 1))
        if self.torsions is not None:
            for atom_1, atom_2, atom_3, atom_4 in self.torsions:
                combined.append("TORSION {} {} {} {}".format(atom_1 + 1,
                                                             atom_2 + 1,
                                                             atom_3 + 1,
                                                             atom_4 + 1))
        if self.out_of_planes is not None:
            for atom_1, atom_2, atom_3, atom_4 in self.out_of_planes:
                combined.append("OOP {} {} {} {}".format(atom_1 + 1,
                                                         atom_2 + 1,
                                                         atom_3 + 1,
                                                         atom_4 + 1))

        return "\n".join(combined)

    def write_file(self, filename):
        """
        Write the isomers file.

        Args:
            filename (str)

        Returns:
            None
        """
        with zopen(filename, 'wt') as f:
            f.write(self.__str__())

    @classmethod
    def from_string(cls, string):
        """
        Extract coordinate information from a string and convert it to a
        GSMIsomerInput object.

        NOTE: This function does no checks to verify that the given string
            follows the appropriate format.

        Args:
            string (str): string to be converted.
        Returns:
            GSMIsomerInput
        """

        lines = string.split("\n")

        bonds_formed = list()
        bonds_broken = list()
        angles = list()
        torsions = list()
        out_of_planes = list()

        for line in lines:
            tokens = line.rstrip().split(" ")
            if tokens[0] == "NEW":
                continue
            elif tokens[0] == "ADD":
                bonds_formed.append((int(tokens[1]) - 1, int(tokens[2]) - 1))
            elif tokens[0] == "BREAK":
                bonds_broken.append((int(tokens[1]) - 1, int(tokens[2]) - 1))
            elif tokens[0] == "ANGLE":
                angles.append((int(tokens[1]) - 1,
                               int(tokens[2]) - 1,
                               int(tokens[3]) - 1))
            elif tokens[0] == "TORSION":
                torsions.append((int(tokens[1]) - 1, int(tokens[2]) - 1,
                                 int(tokens[3]) - 1, int(tokens[4]) - 1))
            elif tokens[0] == "OOP":
                out_of_planes.append((int(tokens[1]) - 1, int(tokens[2]) - 1,
                                      int(tokens[3]) - 1, int(tokens[4]) - 1))

        return cls(bonds_formed=bonds_formed, bonds_broken=bonds_broken,
                   angles=angles, torsions=torsions,
                   out_of_planes=out_of_planes)

    @staticmethod
    def from_file(filename):
        with zopen(filename, 'rt') as f:
            return GSMIsomerInput.from_string(f.read())


def parse_multi_xyz(filename):
    """
    Extract multiple molecules from an XYZ file

    Note: This file will fail if not given a valid XYZ file

    TODO: Do some more elegant parsing to ensure that the xyz file is valid

    Args:
        filename (str): The multi-XYZ file to be parsed.

    Returns:
        molecules (list of Molecule objects)
    """

    molecules = list()

    with open(filename) as molfile:
        text = molfile.readlines()

        linenum = 0

        while linenum < len(text):
            try:
                num_atoms = int(text[linenum].strip())
                mol = Molecule.from_str("".join(text[linenum:linenum + num_atoms + 2]), "xyz")
                molecules.append(mol)

                linenum += num_atoms + 2
            except ValueError:
                break

    return molecules
