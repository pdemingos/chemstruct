# Copyright 2019 Pedro G. Demingos

"""
Chemical Structure Analysis (ChemStruct)
----------------------------------------

Provides
    1. Atom and Atoms objects for atomic representation
    2. Methods for chemical structure analysis of atomic structures
    3. Objects for reading/writing files (e.g. xyz and lmp files)

More tools are currently under development.
For updates, see:
https://github.com/pdemingos/chemstruct

"""

name = "chemstruct"

import numpy as np

# sets this path
from os.path import dirname, abspath
import sys
sys.path.insert(0, dirname(abspath(__file__)))

