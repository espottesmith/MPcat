# coding: utf-8

from typing import Optional
import datetime
from pathlib import Path

from monty.json import MSONable, jsanitize

from pymatgen.core.periodic_table import DummySpecie

from schrodinger.application.jaguar.output import JaguarOutput
from schrodinger.application.jaguar.results import JaguarResults

from mpcat.adapt.schrodinger_adapter import (maestro_file_to_molecule, schrodinger_struct_to_molecule)


class JaguarOutputParseError(Exception):
    """Exception for errors that occur while parsing Jaguar output files"""
    pass


def parse_jaguar_results(jagresult: JaguarResults):
    """
    Helper function to parse a JaguarResults object.

    Args:
        jagresult (JaguarResults object): Results from a Jaguar calculation.

    Returns:
        data (dict): data from Jaguar results
    """

    data = dict()

    data["atoms"] = list()
    for atom in jagresult.atom:
        atom_data = dict()
        atom_data["name"] = atom.atom_name
        atom_data["forces"] = atom.forces
        atom_data["mulliken_charge"] = atom.charge_mulliken
        atom_data["mulliken_spin"] = atom.spin_mulliken
        atom_data["esp_charge"] = atom.charge_esp
        data["atoms"].append(atom_data)

    data["scf_energy"] = jagresult.scf_energy
    data["gas_phase_energy"] = jagresult.gas_phase_energy
    data["solution_phase_energy"] = jagresult.solution_phase_energy
    data["one_electron_energy"] = jagresult.energy_one_electron
    data["two_electron_energy"] = jagresult.energy_two_electron
    data["electronic_energy"] = jagresult.energy_electronic
    data["aposteri_energy"] = jagresult.energy_aposteri
    data["nuclear_repulsion_energy"] = jagresult.nuclear_repulsion
    data["zero_point_energy"] = jagresult.zero_point_energy

    data["homo_alpha"] = jagresult.homo_alpha
    data["homo_beta"] = jagresult.homo_beta
    data["lumo_alpha"] = jagresult.lumo_alpha
    data["lumo_beta"] = jagresult.lumo_beta

    data["reaction_coord"] = jagresult.reaction_coord

    data["thermo"] = list()
    if isinstance(jagresult.thermo, list):
        thermos = jagresult.thermo
    else:
        thermos = [jagresult.thermo]

    for thermo_collection in thermos:
        this_thermo = dict()
        this_thermo["temperature"] = thermo_collection.temp
        this_thermo["pressure"] = thermo_collection.press
        this_thermo["units"] = thermo_collection.units
        this_thermo["total_energy"] = thermo_collection.UTotal
        this_thermo["total_enthalpy"] = thermo_collection.HTotal
        this_thermo["total_free_energy"] = thermo_collection.GTotal

        this_thermo["energy"] = dict()
        this_thermo["energy"]["temperature"] = thermo_collection.U.temp
        this_thermo["energy"]["pressure"] = thermo_collection.U.press
        this_thermo["energy"]["units"] = thermo_collection.U.energy_units
        this_thermo["energy"]["total_energy"] = thermo_collection.U.total
        this_thermo["energy"]["translational_energy"] = thermo_collection.U.translational
        this_thermo["energy"]["rotational_energy"] = thermo_collection.U.rotational
        this_thermo["energy"]["vibrational_energy"] = thermo_collection.U.vibrational
        this_thermo["energy"]["electronic_energy"] = thermo_collection.U.electronic

        this_thermo["enthalpy"] = dict()
        this_thermo["enthalpy"]["temperature"] = thermo_collection.H.temp
        this_thermo["enthalpy"]["pressure"] = thermo_collection.H.press
        this_thermo["enthalpy"]["units"] = thermo_collection.H.energy_units
        this_thermo["enthalpy"]["total_enthalpy"] = thermo_collection.H.total
        this_thermo["enthalpy"]["translational_enthalpy"] = thermo_collection.H.translational
        this_thermo["enthalpy"]["rotational_enthalpy"] = thermo_collection.H.rotational
        this_thermo["enthalpy"]["vibrational_enthalpy"] = thermo_collection.H.vibrational
        this_thermo["enthalpy"]["electronic_enthalpy"] = thermo_collection.H.electronic

        this_thermo["entropy"] = dict()
        this_thermo["entropy"]["temperature"] = thermo_collection.S.temp
        this_thermo["entropy"]["pressure"] = thermo_collection.S.press
        this_thermo["entropy"]["units"] = thermo_collection.S.energy_units
        this_thermo["entropy"]["total_entropy"] = thermo_collection.S.total
        this_thermo["entropy"]["translational_entropy"] = thermo_collection.S.translational
        this_thermo["entropy"]["rotational_entropy"] = thermo_collection.S.rotational
        this_thermo["entropy"]["vibrational_entropy"] = thermo_collection.S.vibrational
        this_thermo["entropy"]["electronic_entropy"] = thermo_collection.S.electronic

        this_thermo["free_energy"] = dict()
        this_thermo["free_energy"]["temperature"] = thermo_collection.G.temp
        this_thermo["free_energy"]["pressure"] = thermo_collection.G.press
        this_thermo["free_energy"]["units"] = thermo_collection.G.energy_units
        this_thermo["free_energy"]["total_free_energy"] = thermo_collection.G.total
        this_thermo["free_energy"]["translational_free_energy"] = thermo_collection.G.translational
        this_thermo["free_energy"]["rotational_free_energy"] = thermo_collection.G.rotational
        this_thermo["free_energy"]["vibrational_free_energy"] = thermo_collection.G.vibrational
        this_thermo["free_energy"]["electronic_free_energy"] = thermo_collection.G.electronic

        this_thermo["c_v"] = dict()
        this_thermo["c_v"]["temperature"] = thermo_collection.Cv.temp
        this_thermo["c_v"]["pressure"] = thermo_collection.Cv.press
        this_thermo["c_v"]["units"] = thermo_collection.Cv.energy_units
        this_thermo["c_v"]["total_heat_capacity"] = thermo_collection.Cv.total
        this_thermo["c_v"]["translational_heat_capacity"] = thermo_collection.Cv.translational
        this_thermo["c_v"]["rotational_heat_capacity"] = thermo_collection.Cv.rotational
        this_thermo["c_v"]["vibrational_heat_capacity"] = thermo_collection.Cv.vibrational
        this_thermo["c_v"]["electronic_heat_capacity"] = thermo_collection.Cv.electronic

        this_thermo["partition_function"] = dict()
        this_thermo["partition_function"]["temperature"] = thermo_collection.lnQ.temp
        this_thermo["partition_function"]["pressure"] = thermo_collection.lnQ.press
        this_thermo["partition_function"]["units"] = thermo_collection.lnQ.energy_units
        this_thermo["partition_function"]["total_partition_function"] = thermo_collection.lnQ.total
        this_thermo["partition_function"]["translational_partition_function"] = thermo_collection.lnQ.translational
        this_thermo["partition_function"]["rotational_partition_function"] = thermo_collection.lnQ.rotational
        this_thermo["partition_function"]["vibrational_partition_function"] = thermo_collection.lnQ.vibrational
        this_thermo["partition_function"]["electronic_partition_function"] = thermo_collection.lnQ.electronic

        data["thermo"].append(this_thermo)

    data["dipole"] = dict()
    data["dipole"]["qm"] = {"magnitude": jagresult.dipole_qm.magnitude,
                            "vector": [jagresult.dipole_qm.x,
                                       jagresult.dipole_qm.y,
                                       jagresult.dipole_qm.z]}
    data["dipole"]["esp"] = {"magnitude": jagresult.dipole_esp.magnitude,
                             "vector": [jagresult.dipole_esp.x,
                                        jagresult.dipole_esp.y,
                                        jagresult.dipole_esp.z]}
    data["dipole"]["mulliken"] = {"magnitude": jagresult.dipole_mulliken.magnitude,
                                  "vector": [jagresult.dipole_mulliken.x,
                                             jagresult.dipole_mulliken.y,
                                             jagresult.dipole_mulliken.z]}

    data["frequencies"] = list()
    data["vibrational_frequency_modes"] = list()
    if isinstance(jagresult.normal_mode, list):
        modes = jagresult.normal_mode
    else:
        modes = [jagresult.normal_mode]

    try:
        data["scan_values"] = jagresult.scan_values
    except AttributeError:
        data["scan_values"] = None

    for mode in modes:
        data["frequencies"].append(mode.frequency)
        data["vibrational_frequency_modes"].append(mode.displacement)

    data["rotational_constants"] = jagresult.rotational_constants

    molecule = schrodinger_struct_to_molecule(jagresult.getStructure())
    data["molecule"] = molecule

    return data


class JagOutput(MSONable):
    """
    Class to parse Jaguar output files and convert from Schrodinger-specific
    data types/objects to standard Python/pymatgen types/objects.
    """

    def __init__(self, filename: str, allow_failure: Optional[bool] = False,
                 parse_molecules: Optional[bool] = True):
        """
        Args:
            filename (str): Path to Jaguar output file to be parsed
            allow_failure (bool): If True (default False), do not raise an
                exception if parsing fails. Capture what information can be
                parsed.
            parse_molecules (bool): If True (default), then in addition to
                parsing the Jaguar output file provided by filename, any
                molecule (*.mae) files referenced in this output file will be
                parsed.

        Returns:
            None
        """

        self.filename = filename
        self.data = dict()
        if self.filename != "":
            base_dir = Path(self.filename).resolve().parent
            try:
                jag_out = JaguarOutput(output=filename, partial_ok=allow_failure)
            except StopIteration:
                raise JaguarOutputParseError

            self.data["job_name"] = jag_out.name
            self.data["full_filename"] = jag_out.filename
            self.data["job_id"] = jag_out.job_id
            self.data["parsing_error"] = jag_out.parsing_error
            self.data["fatal_error"] = jag_out.fatal_error
            self.data["point_group"] = jag_out.point_group

            if jag_out.status == 1:
                self.data["success"] = True
            else:
                self.data["success"] = False

            self.data["start_time"] = jag_out.start_time
            self.data["end_time"] = jag_out.end_time

            if self.data["start_time"] is not None and self.data["end_time"] is not None:
                self.data["walltime"] = (self.data["end_time"] - self.data["start_time"]).total_seconds()
            else:
                self.data["walltime"] = None

            self.data["input"] = dict()
            self.data["input"]["basis"] = jag_out.basis
            self.data["input"]["functional"] = jag_out.functional
            self.data["input"]["method"] = jag_out.method
            self.data["input"]["canonical_orbitals"] = jag_out.canonical_orbitals
            self.data["input"]["charge"] = jag_out.charge
            self.data["input"]["multiplicity"] = jag_out.multiplicity
            self.data["input"]["basis_functions"] = jag_out.nbasis
            self.data["input"]["num_electrons"] = jag_out.nelectron
            self.data["input"]["input_file"] = jag_out.mae_in
            self.data["input"]["solvation"] = jag_out.opts.solvation

            self.data["energy_trajectory"] = list()
            for opt_step in jag_out.geopt_step:
                parsed = parse_jaguar_results(opt_step)
                self.data["energy_trajectory"].append(parsed["scf_energy"])

            self.data["output"] = parse_jaguar_results(jag_out._results)
            self.data["output"]["irc"] = [parse_jaguar_results(j) for j in jag_out.irc_step]
            self.data["output"]["scan"] = [parse_jaguar_results(j) for j in jag_out.scan_step]
            self.data["output"]["output_file"] = jag_out.mae_out

            if parse_molecules:
                self.data["input"]["molecule"] = maestro_file_to_molecule(
                    base_dir / self.data["input"]["input_file"])[0]
                self.data["molecule_trajectory"] = maestro_file_to_molecule(
                    base_dir / self.data["output"]["output_file"])
                self.data["output"]["molecule"] = self.data["molecule_trajectory"][-1]

                self.data["input"]["molecule"].remove_species([DummySpecie("")])
                for m in self.data["molecule_trajectory"]:
                    m.remove_species([DummySpecie("")])
                self.data["output"]["molecule"].remove_species([DummySpecie("")])

            else:
                self.data["input"]["molecule"] = None
                self.data["output"]["molecule"] = None
                self.data["molecule_trajectory"] = None

    def as_dict(self):
        d = dict()
        d["data"] = self.data
        try:
            d["data"]["start_time"] = self.data["start_time"].strftime("%Y/%m/%d, %H:%M:%S")
        except AttributeError:
            d["data"]["start_time"] = None
        try:
            d["data"]["end_time"] = self.data["end_time"].strftime("%Y/%m/%d, %H:%M:%S")
        except AttributeError:
            d["data"]["end_time"] = None

        d["filename"] = self.filename
        return jsanitize(d, strict=True)

    @classmethod
    def from_dict(cls, d):
        output = JagOutput("", allow_failure=True,
                           parse_molecules=False)
        output.data = d["data"]
        if d["data"]["start_time"] is not None:
            output.data["start_time"] = datetime.datetime.strptime(d["data"]["start_time"],
                                                                   "%Y/%m/%d, %H:%M:%S")
        if d["data"]["end_time"] is not None:
            output.data["end_time"] = datetime.datetime.strptime(d["data"]["end_time"],
                                                                 "%Y/%m/%d, %H:%M:%S")
        output.filename = d["filename"]

        return output
