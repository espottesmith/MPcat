from monty.io import zopen
from monty.json import MSONable
from pymatgen.io.qchem.utils import (read_table_pattern,
                                     read_pattern)


class AutoTSOutput(MSONable):
    """
    Class to parse AutoTS output files.
    """

    def __init__(self, filename: str):
        """
        Args:
            filename (str): Path to AutoTS output file.
        """

        self.filename = filename
        self.data = dict()
        self.data["warnings"] = dict()
        self.data["input"] = dict()
        self.data["output"] = dict()
        self.data["errors"] = dict()
        self.text = ""
        with zopen(filename, 'rt') as f:
            self.text = f.read()

        self._parse_input()

        self._parse_calculations()

        self._parse_warnings()
        self._parse_errors()

        calculation_summary = read_pattern(self.text,
                                           {"key": r"Calculation Summary"},
                                           terminate_on_match=True).get("key")

        time = read_pattern(self.text,
                            {"key": r"Timer \(Total AutoTS Time\) : ([0-9\.]+) secs \([0-9A-Za-z,\. ]+\)"},
                            terminate_on_match=True).get("key")
        if time is not None:
            self.data["walltime"] = float(time[0][0])
        else:
            self.data["walltime"] = None

        if calculation_summary is None:
            self.data["complete"] = False
        else:
            self.data["complete"] = True
            self._parse_calculation_summary()

    def _parse_input(self):
        version = read_pattern(
            self.text, {"key": r"\s+Version: ([0-9\-]+)"}, terminate_on_match=True).get('key')
        self.data["input"]["version"] = version[0][0]

        schro = read_pattern(self.text, {"key": r"\s+\$SCHRODINGER: ([A-Za-z0-9_\.\-/]+)"},
                             terminate_on_match=True).get("key")
        self.data["input"]["schrodinger_path"] = schro[0][0]

        launch_dir = read_pattern(self.text, {"key": r"\s+Launch Directory: ([A-Za-z0-9_\.\-/]+)"},
                                  terminate_on_match=True).get("key")
        self.data["input"]["launch_directory"] = launch_dir[0][0]
        run_dir = read_pattern(self.text, {"key": r"\s+Run Directory: ([A-Za-z0-9_\.\-/]+)"},
                               terminate_on_match=True).get("key")
        self.data["input"]["run_directory"] = run_dir[0][0]

        lot = read_pattern(self.text, {"key": r"\s+Level of Theory: ([A-Za-z0-9\-]+)/([A-Za-z0-9\-]+)"},
                           terminate_on_match=True).get("key")
        self.data["input"]["functional"] = lot[0][0]
        self.data["input"]["basis"] = lot[0][1]
        solution = read_pattern(self.text, {"key": r"\sLevel of Theory: .*, Solution Phase"},
                                terminate_on_match=True).get("key")
        if solution is None:
            self.data["input"]["solution_phase"] = False
        else:
            self.data["input"]["solution_phase"] = True

    def _parse_calculations(self):
        temp_calcs = read_pattern(self.text,
                                  {"key": r"\([0-9]+\) jaguar run ([A-Za-z0-9_\.\-]+) -TPP [0-9]+"}).get("key")

        if temp_calcs is None or len(temp_calcs) == 0:
            self.data["calculations"] = None
        else:
            self.data["calculations"] = list()
            for calc in temp_calcs:
                self.data["calculations"].append(calc[0].replace(".in", ""))

    def _parse_warnings(self):
        temp_warnings = read_pattern(self.text,
                                     {"syn_addition": r"Warning: Not performing syn\-addition correction",
                                      "ff_typing": r"Warning: FF typing failed, cannot add VDW in path optimization",
                                      "complexes": r"Warning: Unable to optimize initial complexes",
                                      "canonical": (r"Warning: Exception raised while locating transition states:" +
                                                    r" cannot determine how many canonical orbitals to take"),
                                      "lewis": r"mmlewis warning"},
                                     terminate_on_match=True)
        if temp_warnings.get("syn_addition"):
            self.data["warnings"]["syn_addition"] = True
        if temp_warnings.get("ff_typing"):
            self.data["warnings"]["ff_typing"] = True
        if temp_warnings.get("complexes"):
            self.data["warnings"]["complexes"] = True
        if temp_warnings.get("canonical"):
            self.data["warnings"]["canonical"] = True
        if temp_warnings.get("lewis"):
            self.data["warnings"]["lewis"] = True

    def _parse_errors(self):
        temp_errors = read_pattern(self.text,
                                   {"lewis": r"FATAL mmlewis error: Problems remain with the following atoms:",
                                    "canonical": (r"A fatal error occurred: cannot determine how many canonical " +
                                                  r"orbitals to take")},
                                   terminate_on_match=True)

        if temp_errors.get("lewis"):
            self.data["errors"]["lewis"] = True
        if temp_errors.get("canonical"):
            self.data["errors"]["canonical"] = True

        failed_jobs = read_pattern(self.text,
                                   {"key": r"Job failure\(s\)\nJobname: ([A-Za-z0-9\-\._]+)"})
        if failed_jobs.get("key") is not None:
            self.data["failed_jobs"] = list()
            for failed_job in failed_jobs.get("key"):
                self.data["failed_jobs"].append(failed_job[0])

    def _parse_calculation_summary(self):

        header_struct_energies = (r"\s*=+\n\s+Structure Energies \(Hartree\)\s+\n\s*=+\n\s+Structure\s+Total\s+" +
                                  r"Separated\s+Gibbs \(298.15K\)\s*\n\s*=+")
        table_struct_energies = r"\s*([A-Za-z0-9_]+)\s+([\-\.0-9]+)\s+([\-\.0-9]+)\s+([\-\.0-9]+)"
        footer_struct_energies = r"\s*=+"

        temp_struct_energies = read_table_pattern(self.text, header_struct_energies,
                                                  table_struct_energies,
                                                  footer_struct_energies)
        if temp_struct_energies is None or len(temp_struct_energies) == 0:
            self.data["output"]["struct_energies"] = None
        else:
            self.data["output"]["struct_energies"] = dict()
            for row in temp_struct_energies[0]:
                struct_name = row[0]
                total_energy = float(row[1])
                separated_energy = float(row[2])
                gibbs = float(row[3])
                self.data["output"]["struct_energies"][struct_name] = {"total": total_energy,
                                                                       "separated": separated_energy,
                                                                       "gibbs": gibbs}

        header_complex_energy = (r"\s+Complex Energy \(kcal/mol\)\s+\n\s*\-+\n\s*Rxn\s+Er\s+Ep\s+Ets\s+Rxn E\s+Act E" +
                                 r"\s+frequency\s+Status of IRC\?\s+Path number\s+Approx TS")
        table_complex_energy = (r"\s*([0-9]+)\s+([\-\.0-9]+)\s+([\-\.0-9]+)\s+([\-\.0-9]+)\s+([\-\.0-9]+)\s+" +
                                r"([\-\.0-9]+)\s+([\-\.0-9]+)\s+(Failed|Passed)\s+([0-9]+)\s+(True|False)\s*")
        footer_complex_energy = r"\s*\-+"

        temp_complex_energies = read_table_pattern(self.text, header_complex_energy,
                                                   table_complex_energy,
                                                   footer_complex_energy)

        if temp_complex_energies is None or len(temp_complex_energies) == 0:
            self.data["output"]["complex_energies"] = None
        else:
            self.data["output"]["complex_energies"] = list()
            for row in temp_complex_energies[0]:
                if row[7] == "Failed":
                    irc = False
                else:
                    irc = True

                if row[9] == "False":
                    approx_ts = False
                else:
                    approx_ts = True

                self.data["output"]["complex_energies"].append({"reactant": float(row[1]),
                                                                "product": float(row[2]),
                                                                "transition_state": float(row[3]),
                                                                "rxn": float(row[4]),
                                                                "activation": float(row[5]),
                                                                "freq_mode": float(row[6]),
                                                                "irc_succeeded": irc,
                                                                "path_number": int(row[8]),
                                                                "approx_ts": approx_ts
                                                                })

        header_sep_energy = (r"\s+Inf\. Separated Energy \(kcal/mol\)\s+\n\s*\-+\n\s*Rxn\s+Er\s+Ep\s+Ets\s+Rxn E\s+" +
                             r"Act E\s+frequency\s+Status of IRC\?\s+Path number\s+Approx TS")
        table_sep_energy = (r"\s*([0-9]+)\s+([\-\.0-9]+)\s+([\-\.0-9]+)\s+([\-\.0-9]+)\s+([\-\.0-9]+)\s+" +
                            r"([\-\.0-9]+)\s+([\-\.0-9]+)\s+(Failed|Passed)\s+([0-9]+)\s+(True|False)\s*")
        footer_sep_energy = r"\s*\-+"

        temp_sep_energies = read_table_pattern(self.text, header_sep_energy,
                                               table_sep_energy,
                                               footer_sep_energy)
        if temp_sep_energies is None or len(temp_sep_energies) == 0:
            self.data["output"]["separated_energies"] = None
        else:
            self.data["output"]["separated_energies"] = list()
            for row in temp_sep_energies[0]:
                if row[7] == "Failed":
                    irc = False
                else:
                    irc = True

                if row[9] == "False":
                    approx_ts = False
                else:
                    approx_ts = True

                self.data["output"]["separated_energies"].append({"reactant": float(row[1]),
                                                                  "product": float(row[2]),
                                                                  "transition_state": float(row[3]),
                                                                  "rxn": float(row[4]),
                                                                  "activation": float(row[5]),
                                                                  "freq_mode": float(row[6]),
                                                                  "irc_succeeded": irc,
                                                                  "path_number": int(row[8]),
                                                                  "approx_ts": approx_ts
                                                                  })

        header_gibbs_energy = (r"\s+Gibbs Free Energy at 298\.15K and concentration of 1 mol/liter \(kcal/mol\)\s+\n" +
                               r"\s*\-+\n\s*Rxn\s+Gr\s+Gp\s+Gts\s+Rxn G\s+Act G" +
                               r"\s+frequency\s+Status of IRC\?\s+Path number\s+Approx TS")
        table_gibbs_energy = (r"\s*([0-9]+)\s+([\-\.0-9]+)\s+([\-\.0-9]+)\s+([\-\.0-9]+)\s+([\-\.0-9]+)\s+" +
                              r"([\-\.0-9]+)\s+([\-\.0-9]+)\s+(Failed|Passed)\s+([0-9]+)\s+(True|False)\s*")
        footer_gibbs_energy = r"\s*\-+"

        temp_gibbs_energies = read_table_pattern(self.text, header_gibbs_energy,
                                                 table_gibbs_energy,
                                                 footer_gibbs_energy)
        if temp_gibbs_energies is None or len(temp_gibbs_energies) == 0:
            self.data["output"]["gibbs_energies"] = None
        else:
            self.data["output"]["gibbs_energies"] = list()
            for row in temp_gibbs_energies[0]:
                if row[7] == "Failed":
                    irc = False
                else:
                    irc = True

                if row[9] == "False":
                    approx_ts = False
                else:
                    approx_ts = True

                self.data["output"]["gibbs_energies"].append({"reactant": float(row[1]),
                                                                "product": float(row[2]),
                                                                "transition_state": float(row[3]),
                                                                "rxn": float(row[4]),
                                                                "activation": float(row[5]),
                                                                "freq_mode": float(row[6]),
                                                                "irc_succeeded": irc,
                                                                "path_number": int(row[8]),
                                                                "approx_ts": approx_ts
                                                                })

        header_temp_dep = (r"\(Celsius\)\s+\(kcal/mol\)\s+\(kcal/mol\)\s+\(s\^\-1.*\)\s+\(s\^\-1.*\)\s+\(unitless\)" +
                           r"\s+at t1/2 \(%\)\s+half life\s+\n\-+")
        table_temp_dep = (r"\s*([0-9\.\-]+)\s+([0-9\+\-eE\.]+)\s+([0-9\+\-eE\.]+)\s+([0-9\+\-eE\.]+)\s+" +
                          r"([0-9\+\-eE\.]+)\s+([0-9\+\-eE\.]+)\s+([0-9\+\-eE\.]+)\s+([A-Za-z0-9\.\-\+ ]+)")
        footer_temp_dep = r""

        temp_temp_dep = read_table_pattern(self.text, header_temp_dep,
                                           table_temp_dep, footer_temp_dep)
        if temp_temp_dep is None or len(temp_temp_dep) == 0:
            self.data["output"]["temperature_dependence"] = None
        else:
            self.data["output"]["temperature_dependence"] = dict()
            for row in temp_temp_dep[0]:
                temp = float(row[0])
                self.data["output"]["temperature_dependence"][temp] = {
                    "gibbs_activation": float(row[1]),
                    "gibbs_reaction": float(row[2]),
                    "k_forwards": float(row[3]),
                    "k_reverse": float(row[4]),
                    "log10_keq": float(row[5]),
                    "rct_remaining_t_half": float(row[6]),
                    "half_life": row[7]
                }

        header_diagnostic = r"\s+AutoTS diagnostic report\s+\n\-+\n\s+property\s+value\s+diagnostic passes\n\-+"
        table_diagnostic = r"\s*((?:[A-za-z]+ ?)+)\s*([A-Za-z0-9\-,/ ]+)\s+(Yes|No)"
        footer_diagnostic = r""

        temp_diagnostic = read_table_pattern(self.text, header_diagnostic,
                                             table_diagnostic,
                                             footer_diagnostic)
        if temp_diagnostic is None or len(temp_diagnostic) == 0:
            self.data["diagnostic"] = None
        else:
            self.data["diagnostic"] = dict()
            for row in temp_diagnostic[0]:
                if "Level of theory" in row[0]:
                    continue
                if row[2].strip() == "Yes":
                    passed = True
                else:
                    passed = False

                self.data["diagnostic"][row[0].strip()] = {"value": row[1].strip(),
                                                           "passed": passed}
