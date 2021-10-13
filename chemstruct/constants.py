# Copyright 2019 Pedro G. Demingos

"""Deals with some constants and other informations."""

from os.path import dirname, abspath
PACKAGE_DIR = dirname(abspath(__file__))

BOHR_TO_ANGSTROM = 0.529177249
RY_TO_EV = 13.605698066
EV_PER_ANGSTROM_TO_NANONEWT = 1.60217662

BOND_LENGTHS = dict()
ATOMIC_MASSES = dict()


def axis_to_dim(axis):
    if axis in [0, 1, 2, "0", "1", "2"]:
        return int(axis)
    if not isinstance(axis, str):
        raise KeyError("bad axis: it should be a string")
    axis = axis.lower()
    if axis == "x":
        return 0
    elif axis == "y":
        return 1
    elif axis == "z":
        return 2
    else:
        raise KeyError("bad axis, should be 'x', 'y' or 'z'")


def get_params():
    get_bond_lengths()
    get_atomic_masses()
    # ...


def get_bond_lengths():
    global BOND_LENGTHS
    with open(PACKAGE_DIR + "/bond_lengths.txt", "r") as F:
        line = F.readline()
        while line:
            try:
                atom1, atom2, length = tuple(line.split())
            except ValueError:
                raise TypeError("file bond_lengths.txt has wrong format, "
                                "keep it as 'H H 0.8' in every line")
            else:
                bond12 = atom1 + ":" + atom2
                bond21 = atom2 + ":" + atom1
                try:
                    BOND_LENGTHS[bond12] = float(length)
                    BOND_LENGTHS[bond21] = float(length)
                except TypeError:
                    raise TypeError("bad bond length given in bond_lengths.txt for {} bond, "
                                    "expected a number but got {}".format(bond12, length))
            line = F.readline()


def get_atomic_masses():
    global ATOMIC_MASSES
    with open(PACKAGE_DIR + "/masses.txt", "r") as F:
        line = F.readline()
        while line:
            try:
                atom, mass = tuple(line.split())
            except ValueError:
                raise TypeError("file masses.txt has wrong format, "
                                "keep it as 'C 12.00' in every line")
            else:
                try:
                    ATOMIC_MASSES[atom] = float(mass)
                except TypeError:
                    raise TypeError("bad mass given in masses.txt for {} atom, "
                                    "expected a number but got {}".format(atom, mass))
            line = F.readline()


if __name__ == "__main__":
    pass
else:  # if imported
    get_params()
