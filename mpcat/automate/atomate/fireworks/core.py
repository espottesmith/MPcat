# coding: utf-8

# Defines standardized Fireworks that can be chained easily to perform
# calculations using pyGSM

from mpcat.automate.atomate.sets import GSMDictSet

from fireworks import Firework

from mpcat.automate.atomate.firetasks.write_inputs import WriteInputFromIOSet, WriteIsomer
from mpcat.automate.atomate.firetasks.run_calc import RunGSMDirect
from mpcat.automate.atomate.firetasks.parse_outputs import GSMToDb

__author__ = "Evan Spotte-Smith"
__copyright__ = "Copyright 2021"
__version__ = "0.1"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "espottesmith@gmail.com"
__status__ = "Alpha"
__date__ = "4/24/20"
__credits__ = "Sam Blau"


class SingleEndedGSMFW(Firework):
    def __init__(self,
                 molecule,
                 isomers,
                 base_command="gsm",
                 num_nodes=30,
                 fixed_endpoints=True,
                 name="single-ended GSM",
                 max_cores=32,
                 input_params=None,
                 db_file=None,
                 parents=None,
                 **kwargs):
        """

        Args:
            molecule (Molecule): Input molecule.
            isomers (dict): list of coordinates to vary to find transition states
            base_command (str): Way to call pyGSM. By default this is "gsm".
            num_nodes (int): Maximum number of nodes along the reaction path.
                 Default is 30.
            fixed_endpoint (bool): If True (default), do not optimize the input
                node.
            name (str): Name for the Firework.
            max_cores (int): Maximum number of cores to parallelize over. Supports env_chk.
            input_params (dict): Specify kwargs for instantiating the input set parameters.
                                 Basic uses would be to modify the default inputs of the set,
                                 such as dft_rung, basis_set, pcm_dielectric, scf_algorithm,
                                 or max_scf_cycles. See pymatgen/io/gsm/sets.py for default
                                 values of all input parameters. For instance, if a user wanted
                                 to use a more advanced DFT functional, include a pcm with a
                                 dielectric of 30, and use a larger basis, the user would set
                                 qchem_input_params = {"dft_rung": 5, "pcm_dielectric": 30,
                                 "basis_set": "6-311++g**"}. However, more advanced customization
                                 of the input is also possible through the overwrite_inputs key
                                 which allows the user to directly modify the rem, pcm, smd, and
                                 solvent dictionaries that QChemDictSet passes to inputs.py to
                                 print an actual input file. For instance, if a user wanted to
                                 set the sym_ignore flag in the rem section of the input file
                                 to true, then they would set qchem_input_params = {"overwrite_inputs":
                                 "rem": {"sym_ignore": "true"}}. Of course, overwrite_inputs
                                 could be used in conjuction with more typical modifications.
            db_file (str): Path to file specifying db credentials to place output parsing.
            parents ([Firework]): Parents of this particular Firework.
            **kwargs: Other kwargs that are passed to Firework.__init__.
        """

        input_params = input_params or dict()

        molecule_file = "input.xyz"
        lot_file = "qin"
        isomers_file = "isomers.txt"
        output_file = "gsm.out"

        flags = [base_command, "-mode", "SE_GSM", "-xyzfile", molecule_file, "-package",
                 "QChem", "-lot_inp_file", lot_file, "-num_nodes", str(num_nodes),
                 "-isomers", isomers_file, "-charge", str(molecule.charge),
                 "-multiplicity", str(molecule.spin_multiplicity), "-nproc", str(max_cores)]
        if fixed_endpoints:
            flags.append("-reactant_geom_fixed")
        flags.append(">")
        flags.append(output_file)

        t = list()

        input_set = GSMDictSet(**input_params)

        t.append(
            WriteInputFromIOSet(
                input_set=input_set,
                molecules=[molecule],
                molecule_file=molecule_file,
                lot_file=lot_file))
        t.append(
            WriteIsomer(
                isomers=isomers,
                isomers_file=isomers_file))
        t.append(
            RunGSMDirect(
                cmd=" ".join(flags)))
        t.append(
            GSMToDb(
                db_file=db_file,
                molecule_file=molecule_file,
                template_file=lot_file,
                output_file=output_file,
                isomers_file=isomers_file,
                additional_fields={"task_label": name,
                                   "true_charge": molecule.charge,
                                   "special_run_type": "single_ended_gsm"}))
        super(SingleEndedGSMFW, self).__init__(
            t,
            parents=parents,
            name=name,
            **kwargs)