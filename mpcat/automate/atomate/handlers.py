# coding: utf-8

from __future__ import unicode_literals, division

# This module implements new error handlers for GSM runs.

import os

from pymatgen.io.gsm.outputs import GSMOutput
from pymatgen.io.qchem.outputs import QCOutput

from custodian.custodian import ErrorHandler
from custodian.utils import backup


__author__ = "Evan Spotte-Smith"
__copyright__ = "Copyright 2020, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "espottesmith@gmail.com"
__status__ = "Alpha"
__date__ = "02/17/20"


class GSMErrorHandler(ErrorHandler):
    """
    Master error handler for pyGSM runs.
    """

    is_monitor = False

    def __init__(self,
                 command,
                 xyz_file="input.xyz",
                 input_file="gsm.inp",
                 output_file="gsm.out",
                 isomers_file="isomers.txt"):
        """

        Args:
            command (str): GSM command that was run for this process
            xyz_file (str): File path to input geometry
            input_file (str): File path to level of theory template file
            output_file (str): File path to GSM output file
            isomers_file (str): File path to isomers file. This is only used if
                this calculation uses one of the single-ended methods
        """

        self.command = command
        self.args = self.command.split(" ")
        self.xyz_file = xyz_file
        self.input_file = input_file
        self.output_file = output_file
        self.isomers_file = isomers_file

        self.outdata = None
        self.errors = list()
        self.warnings = list()

    def check(self):
        # Check pyGSM output for errors
        self.outdata = GSMOutput(self.output_file).data
        self.errors = self.outdata.get("errors", list())
        self.warnings = self.outdata.get("warnings", list())

        return len(self.errors) > 0

    def correct(self):
        pass

