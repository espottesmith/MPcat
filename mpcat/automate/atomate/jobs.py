# coding: utf-8

from __future__ import unicode_literals, division
import os
import shutil
import subprocess

from custodian.custodian import Job


__author__ = "Evan Spotte-Smith"
__copyright__ = "Copyright 2020, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "espottesmith@gmail.com"
__status__ = "Alpha"
__date__ = "02/17/20"


class GSMJob(Job):
    """
    A basic job using the Growing String Method as implemented in pyGSM.
    """

    def __init__(self,
                 command,
                 charge,
                 multiplicity,
                 xyz_file="input.xyz",
                 input_file="gsm.inp",
                 output_file="gsm.out",
                 isomers_file="isomers.txt",
                 mode="DE_GSM",
                 suffix="",
                 string_id=0,
                 calc_loc=None,
                 save_scratch=False,
                 multimode=True,
                 max_cores=32,
                 ends_fixed=False,
                 num_nodes=9,
                 additional_options=None):
        """
        Args:
            command (str): Command to run pyGSM.
            charge (int): Molecule charge for this calculation.
            multiplicity (int): Spin multiplicity for this calculation.
            xyz_file (str): Name of the input geometry file for GSM.
            input_file (str): Name of the electronic structure package input file
                for the calculation (format depends on the package).
            output_file (str): Name of the file where all output from pyGSM
                should be piped.
            isomers_file (str): Name of the isomers file to define which
                coordinates should be adjusted. This is necessary only for the
                single-ended methods like SE_GSM.
            mode (str): One of DE_GSM (for the double-ended growing-string
                method), SE_GSM (for the single-ended growing-string method), or
                SE_Cross (for the crossing single-ended growing-string method).
                Default is "DE_GSM".
            suffix (str): String to append to the file in postprocess.
            string_id (int): ID for scratch file directory
            calc_loc (str): Path where pyGSM should run. Defaults to None, in
                which case the program will run in the system-defined SCRATCH.
            save_scratch (bool): Whether to save full scratch directory contents.
                Defaults to False.
            multimode (bool): If True (default), then use multiprocessing.
            max_cores (int): Number of processes to use for parallelization
            ends_fixed (bool): If True (default False), do not optimize the ends
                of the string (the reactants and products).
            num_nodes (int): Number of nodes along the string. Default is 9 (for
                double-ended string, this is generally appropriate), but should
                be changed to a larger number (perhaps 20-30) for the
                single-ended string method.
            additional_options (dict, or None): If not None (default), this
                contains key-value pairs for all additional command-line options
                for pyGSM.
        """
        # For now, only allow Q-Chem. Eventually can allow others like XTB
        self.package = "QChem"

        self.command = command.split(" ")

        if mode in ["DE_GSM", "SE_GSM", "SE_Cross"]:
            self.mode = mode
        else:
            raise ValueError("Invalid GSM mode given. "
                             "Options include DE_GSM, SE_GSM, and SE_Cross.")
        self.charge = charge
        self.spin_multiplicity = multiplicity
        self.xyz_file = xyz_file
        self.input_file = input_file
        self.output_file = output_file
        self.isomers_file = isomers_file
        self.suffix = suffix
        self.string_id = string_id
        self.calc_loc = calc_loc
        self.save_scratch = save_scratch
        self.multimode = multimode
        self.max_cores = max_cores
        self.ends_fixed = ends_fixed
        self.num_nodes = num_nodes

        if additional_options is None:
            self.additional_options = dict()
        else:
            self.additional_options = additional_options

    @property
    def current_command(self):
        command = self.command + ["-mode", self.mode, "-xyzfile", self.xyz_file,
                                  "-lot_inp_file", self.input_file, "-package",
                                  self.package, "-ID", str(self.string_id),
                                  "-num_nodes", str(self.num_nodes), "-charge",
                                  str(self.charge), "-multiplicity",
                                  self.spin_multiplicity]

        if self.mode in ["SE_GSM", "SE_Cross"]:
            command += ["-isomers", self.isomers_file]

        if self.ends_fixed:
            if self.mode == "DE_GSM":
                com = ["-reactant_geom_fixed", "-product_geom_fixed"]
            else:
                com = ["-reactant_geom_fixed"]
            command += com

        if self.multimode:
            com = ["-nproc", str(self.max_cores)]
            command += com

        for key, value in self.additional_options.items():
            command.append("-" + key)
            if value is not None:
                command.append(value)

        com_str = " ".join(command)

        return com_str

    def setup(self):
        if self.multimode:
            os.environ['OMP_NUM_THREADS'] = str(self.max_cores)
        os.environ["SCRATCH"] = os.getcwd()
        print("Current SCRATCH is: {}".format(os.environ["SCRATCH"]))
        if self.calc_loc is not None:
            os.environ["SCRATCH"] = self.calc_loc

    def postprocess(self):
        scratch_dir = os.path.join(os.environ["SCRATCH"], "scratch")
        if self.suffix != "":
            shutil.move(self.input_file, self.input_file + self.suffix)
            shutil.move(self.output_file, self.output_file + self.suffix)
        if not self.save_scratch:
            for file in os.listdir(scratch_dir):
                if file.endswith(".0") or "rem" in file or "tmp" in file or "zmat" in file:
                    os.remove(os.path.join(scratch_dir, file))

    def run(self):
        """
        Perform the actual GSM run.

        Returns:
            (subprocess.Popen) Used for monitoring.
        """
        local_scratch = os.path.join(os.environ["SCRATCH"], "scratch")
        if os.path.exists(local_scratch):
            shutil.rmtree(local_scratch)
        p = subprocess.Popen(self.current_command,
                             stdout=open(self.output_file, 'w'),
                             shell=True)
        return p
