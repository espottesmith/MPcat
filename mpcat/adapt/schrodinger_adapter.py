# coding: utf-8

import random
import os

from pymatgen.core.structure import Molecule

from schrodinger.structure import Structure, StructureReader, StructureWriter


def file_to_schrodinger_structure(filename: str):
    """
    Convert a file (sdf, pdb, sd, mol2, maestro, or maestro_text) to a
        Schrodinger Structure object.

    Args:
        filename (str): Name of the file

    Returns:
        structures (list of schrodinger.structure.Structure objects)
    """

    reader = StructureReader(filename)

    structures = [s for s in reader]
    return structures


def maestro_file_to_molecule(filename: str):
    """
    Convert a Maestro file (*.mae) to pymatgen Molecule objects

    Args:
        filename (str): Name of the file to be converted

    Returns:
        molecules (list of pymatgen.core.structure.Molecule objects)
    """
    structures = file_to_schrodinger_structure(filename=filename)

    molecules = [schrodinger_struct_to_molecule(s) for s in structures]
    return molecules


def molecule_to_maestro_file(molecule: Molecule, filename: str):
    """
    Write a Pymatgen Molecule object to a Maestro file.

    Args;
        molecule (Molecule): Molecule to be written
        filename (str): Path to file where molecule will be written.

    Returns:
         None
    """

    struct = molecule_to_schrodinger_struct(molecule)
    StructureWriter.write(struct, filename)


def schrodinger_struct_to_molecule(structure: Structure):
    """
    Convert a Structure object from Schrodinger to a pymatgen Molecule object.

    Args:
        structure (schrodinger.structure.Structure object): Structure to be
            converted

    Returns:
        mol: pymatgen.core.structure.Molecule object
    """

    formal_charge = structure.formal_charge

    elements = list()
    positions = list()
    for molecule in structure.molecule:
        for atom in molecule.atom:
            elements.append(atom.element)
            positions.append(atom.xyz)

    mol = Molecule(elements, positions)
    mol.set_charge_and_spin(charge=formal_charge)

    return mol


def molecule_to_schrodinger_struct(molecule: Molecule):
    """
    Convert a pymatgen Molecule object to a Schrodinger Structure object

    Args:
        molecule (pymatgen.core.structure.Molecule): Molecule to be converted

    Returns:
        struct: schrodinger.structure.Structure object
    """

    # First need to convert molecule to file to use StructureReader

    file_suffix = random.randint(1, 100000000)
    molecule.to('sdf', "temp_conversion{}.sdf".format(file_suffix))

    reader = StructureReader("temp_conversion{}.sdf".format(file_suffix))

    # Assume only one structure (should be the case for a single SDF file)
    struct = [r for r in reader][0]

    for atom in struct.atom:
        atom.formal_charge = 0
    struct.atom[1].formal_charge = molecule.charge

    os.remove("temp_conversion{}.sdf".format(file_suffix))

    return struct
