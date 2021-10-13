"""Example 01: Writing LAMMPS Data File with atom type 'full'"""

from quick import write_lmpdat

# path to some xyz file
from constants import PACKAGE_DIR  # just the package directory
xyz = PACKAGE_DIR + "/examples/01_alkane.xyz"

# path where the LAMMPS Data File will be writen
lmp = PACKAGE_DIR + "/examples/01_alkane.lmp"

write_lmpdat(xyz, lmp)

# a par file will be writen, but will contain a few None parameters
#
# the function will wait while you complete the par file
# which was writen in the same directory as the input xyz file
#
# here's what you'll have to do:
#
# first, the charges
# (CGenFF has no one-to-one atomtype-charge correspondence, that's why)
# "0.00" is ok for "CG301" since all neighboring groups are neutral
# "-0.09" is ok for "CG311" since its bonded to one H with charge "+0.09"
#
# then, angles and dihedrals that aren't in the CGenFF parameters file:
# "58.35 113.5 11.16 2.561" for "CG331:CG311:CG301"
# since these are the parameters of other C:C:C angles in the par file
# "0.2 3 0 1.0" for all the empty C:C:C:C dihedrals, for the same reason
#
# then, press anything, and the LAMMPS Data File will be writen

# if you press the key before having the par file right, see example 02
