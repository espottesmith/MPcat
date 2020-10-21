
# coding: utf-8

import os
import subprocess
import shutil
from typing import Optional, List, Dict, Type
import datetime

from pymatgen.core.structure import Molecule

from mpcat.apprehend.autots_input import AutoTSSet
from mpcat.utils.comparison import compositions_equal