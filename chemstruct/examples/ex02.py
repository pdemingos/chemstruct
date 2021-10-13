"""Example 02: Finishing LAMMPS Data File with atom type 'full'"""

# first, run example 01 without caring for instructions
# (it runs the function write_lmpdat)

# now press any key (asked by ex01) without editing the par file
# this will give you a LAMMPS Data File full of None parameters
# which of course cannot be used by LAMMPS
#
# to fix that, you can either run write_lmpdat again or do the following

# example 01 will run write_lmpdat and write a par file
# the par file requires editing before being used in the LAMMPS Data File
#
# first, check the charges; some will be None
# (CGenFF has no one-to-one atomtype-charge correspondence, that why)
# use "0.00" for "CG301" since all neighboring groups are neutral
# use "-0.09" for "CG311" since its bonded to one H with charge "+0.09"
#
# then, angles and dihedrals that aren't in the CGenFF parameters file:
# "58.35 113.5 11.16 2.561" for "CG331:CG311:CG301",
# since these are the parameters of other C:C:C angles in the par file
# "0.2 3 0 1.0" for all the empty dihedrals, for the same reason

from quick import finish_lmpdat
from constants import PACKAGE_DIR  # just the package directory

# this is the par file writen previously and already edited
par = PACKAGE_DIR + "/examples/01_alkane(computed_parameters).par"

# path where the LAMMPS Data File will be writen
lmp = PACKAGE_DIR + "/examples/01_alkane.lmp"

# this is the xyz file writen previously, with CGenFF classifications
xyz_types = PACKAGE_DIR + "/examples/01_alkane.types.par"

finish_lmpdat(xyz_types, lmp, par)

# this will give you a good LAMMPS Data File
