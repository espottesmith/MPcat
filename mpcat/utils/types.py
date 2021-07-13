from enum import Enum


class JaguarJobType(Enum):
    SP = 0
    OPT = 1
    TS = 2
    FREQ = 3
    SCAN = 4
    IRC = 5

job_type_mapping = {"sp": JaguarJobType.SP,
                    "opt": JaguarJobType.OPT,
                    "ts": JaguarJobType.TS,
                    "freq": JaguarJobType.FREQ,
                    "scan": JaguarJobType.SCAN,
                    "irc": JaguarJobType.IRC}