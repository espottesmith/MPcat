# coding: utf-8

# This module defines tasks that support running pyGSM (with Q-Chem) in various ways.

import shutil
import os
import subprocess

from pymatgen.core.structure import Molecule
from mpcat.automate.atomate.inputs import QCTemplate, GSMIsomerInput

#TODO: Add custodian support with error handlers

# from custodian import Custodian
# from custodian.qchem.handlers import QChemErrorHandler
# from custodian.gsm.handlers import GSMErrorHandler
# from custodian.qchem.jobs import QCJob

from fireworks import explicit_serialize, FiretaskBase

from atomate.utils.utils import env_chk, get_logger
import numpy as np

__author__ = "Evan Spotte-Smith"
__copyright__ = "Copyright 2021"
__version__ = "0.1"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "espottesmith@gmail.com"
__status__ = "Alpha"
__date__ = "04/17/2020"
__credits__ = "Shyam Dwaraknath, Sam Blau"

logger = get_logger(__name__)


@explicit_serialize
class RunGSMDirect(FiretaskBase):
    """
    Execute a command directly (no custodian).

    Required params:
        cmd (str): The name of the full command line call to run. This should include any
                    for parallelization, saving scratch, and input / output files.
                    Does NOT support env_chk.
    Optional params:
        scratch_dir (str): Path to the scratch directory. Defaults to "/dev/shm/qcscratch/".
                           Supports env_chk.

    """

    required_params = ["cmd"]
    optional_params = ["scratch_dir"]

    def run_task(self, fw_spec):
        cmd = self.get("cmd")
        scratch_dir = env_chk(self.get("scratch_dir"), fw_spec)
        if scratch_dir is None:
            scratch_dir = "/dev/shm/qcscratch/"
        os.putenv("QCSCRATCH", scratch_dir)

        logger.info("Running command: {}".format(cmd))
        return_code = subprocess.call(cmd, shell=True)
        logger.info("Command {} finished running with return code: {}".format(
            cmd, return_code))


@explicit_serialize
class RunNoGSM(FiretaskBase):
    """
    Do NOT run GSM. Do nothing.
    """

    def run_task(self, fw_spec):
        pass


@explicit_serialize
class RunGSMFake(FiretaskBase):
    """
     GSM emulator

     Required params:
         ref_dir (string): Path to reference gsm run directory with appropriate
            inputs and outputs.

     """
    required_params = ["ref_dir"]
    optional_params = ["lot_file", "molecule_file", "isomers_file"]

    def run_task(self, fw_spec):
        self._verify_inputs()
        self._clear_inputs()
        self._generate_outputs()

    def _verify_inputs(self):
        lot_file = self.get("lot_file", "qin")
        user_qin = QCTemplate.from_file(os.path.join(os.getcwd(), "mol.qin"))

        # Check mol.qin
        ref_qin = QCTemplate.from_file(os.path.join(self["ref_dir"], lot_file))

        for key in ref_qin.rem:
            if user_qin.rem.get(key) != ref_qin.rem.get(key):
                raise ValueError("Rem key {} is inconsistent!".format(key))
        if ref_qin.pcm is not None:
            for key in ref_qin.pcm:
                if user_qin.pcm.get(key) != ref_qin.pcm.get(key):
                    raise ValueError("PCM key {} is inconsistent!".format(key))
        if ref_qin.solvent is not None:
            for key in ref_qin.solvent:
                if user_qin.solvent.get(key) != ref_qin.solvent.get(key):
                    raise ValueError(
                        "Solvent key {} is inconsistent!".format(key))
        if ref_qin.smx is not None:
            for key in ref_qin.smx:
                if user_qin.smx.get(key) != ref_qin.smx.get(key):
                    raise ValueError(
                        "SMX key {} is inconsistent!".format(key))

        # Check molecule file
        if self.get("molecule_file"):
            mol_file = self.get("molecule_file")
            user_mol = Molecule.from_file(os.path.join(os.getcwd(),
                                                       "input.xyz"))

            ref_mol = Molecule.from_file(os.path.join(self["ref_dir"],
                                                      mol_file))

            np.testing.assert_equal(ref_mol.species,
                                    user_mol.species)
            np.testing.assert_allclose(
                ref_mol.cart_coords,
                user_mol.cart_coords,
                atol=0.0001)

        if self.get("isomers_file"):
            isomers_file = self.get("isomers_file")
            user_iso = GSMIsomerInput.from_file(os.path.join(os.getcwd(),
                                                             "isomers.txt"))

            ref_iso = GSMIsomerInput.from_file(os.path.join(os.getcwd(),
                                                            isomers_file))

            np.testing.assert_equal(user_iso.num_coords,
                                    ref_iso.num_coords)
            np.testing.assert_equal(user_iso.bonds_formed,
                                    ref_iso.bonds_formed)
            np.testing.assert_equal(user_iso.bonds_broken,
                                    ref_iso.bonds_broken)
            np.testing.assert_equal(user_iso.angles,
                                    ref_iso.angles)
            np.testing.assert_equal(user_iso.torsions,
                                    ref_iso.torsions)
            np.testing.assert_equal(user_iso.out_of_planes,
                                    ref_iso.out_of_planes)

        logger.info("RunQChemFake: verified input successfully")

    @staticmethod
    def _clear_inputs():
        q = os.path.join(os.getcwd(), "qin")
        i = os.path.join(os.getcwd(), "isomers.txt")
        m = os.path.join(os.getcwd(), "input.xyz")
        if os.path.exists(q):
            os.remove(q)
        if os.path.exists(i):
            os.remove(i)
        if os.path.exists(m):
            os.remove(m)

    def _generate_outputs(self):
        # pretend to have run pyGSM by copying pre-generated output from reference dir to cur dir
        for file_name in os.listdir(self["ref_dir"]):
            full_file_name = os.path.join(self["ref_dir"], file_name)
            if os.path.isfile(full_file_name):
                shutil.copy(full_file_name, os.getcwd())
        shutil.copytree(os.path.join(self["ref_dir"], "scratch"),
                        os.path.join(os.getcwd(), "scratch"))
        logger.info("RunQChemFake: ran fake QChem, generated outputs")
