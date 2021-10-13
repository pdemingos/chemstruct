"""Example 03: Writing LAMMPS Data File for a periodic system"""

from quick import write_lmpdat

# path to the xyz file
from constants import PACKAGE_DIR  # just the package directory
xyz = PACKAGE_DIR + "/examples/rubber.xyz"

# you can open rubber.xyz as a text file
# notice how the Lattice is defined in the second line
# the software Ovito can write xyz files with Lattice for you
#
# you can also open it in Ovito or another visualization software
# it's a little chaotic, since it's an amorphous polymer
# (about natural rubber: https://en.wikipedia.org/wiki/Natural_rubber)

# path where the LAMMPS Data File will be writen
lmp = PACKAGE_DIR + "/examples/rubber.lmp"

write_lmpdat(xyz, lmp, periodic="xyz")

# notice that now we specify that the system is periodic in all directions
#
# this will take a while, since there are 2608 atoms in the system
#
# complete the par file with reasonable parameters (see ex01 for examples)

# pressing any key will give you the LAMMPS Data File
