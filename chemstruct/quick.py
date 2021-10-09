# Copyright 2019 Pedro G. Demingos

"""Quick functions for generating LAMMPS Data files."""

from constants import PACKAGE_DIR
from files.ff import CharmmGeneral, GeneralAmber, Compass
from files.lmp import LmpDat


def write_lmpdat_charmm(xyz_file: str, lmpdat_file: str, par_file=None,
                        periodic="", simple=False):
    """
    Writes a LAMMPS Data file with atom_style='full' from xyz file,
    using parameters from the CHARMM General Force Field.
    Requires user to manually check a generated parameters (par) file.

    Parameters
    ----------
    xyz_file : str
        Path to the input xyz file (should include Lattice, as Ovito does).
    lmpdat_file : str
        Path to the output LAMMPS Data file.
    par_file : str, optional
        Path to the input parameters file (.par). If none if given, a file
        will be generated.
    periodic : str, optional
        Axes in which the system is periodic (e.g. 'xyz', 'xy', 'z';
        use an empty string '' for non-periodic). Standard is ''.
    simple : bool, optional
        If a simple, less optimised algorithm is wanted for computing
        bonds, angles and dihedrals. (For small systems. See Notes.)
        Standard is False.

    Notes
    -----
    Writes a par file and asks the user to modify it before pressing any key,
    then finishes the LAMMPS Data File with the parameters. The par file must
    be modified with the parameters not specified by the potential, like
    charges for CGenFF and possibly other structures.

    If a system is periodic in all directions, but one of its dimensions is
    such that it's only divided into 1 or 2 regions, simple=True should be
    used to avoid bugs.

    Examples
    --------
    See directory 'examples'.

    """

    periodic = periodic.lower()

    # CGenFF parameters from the original paper
    # prm_file = PACKAGE_DIR + "/par_all36_cgenff.prm"

    charmm = CharmmGeneral(xyz_file=xyz_file)
    charmm.compute_topology(periodic=periodic, simple=simple)
    if par_file is None:
        par_file = xyz_file.replace(".xyz", "(computed_parameters).par")
    charmm.write_par(xyz_file.replace(".xyz", "(computed_parameters).par"))
    charmm.atoms.write_xyz(xyz_file.replace(".xyz", ".types.xyz"),
                           with_classification=True)

    input("Press anything when the par file is ready.")

    lmp = LmpDat()
    lmp.atoms = charmm.atoms
    lmp.get_params(par_file)
    lmp.cell_to_lo_hi()
    lmp.write_lmpdat(lmpdat_file)


def write_lmpdat_amber(xyz_file: str, lmpdat_file: str, par_file=None,
                       periodic="", simple=False):
    """
    Mostly done. Only improper dihedrals missing.
    """

    periodic = periodic.lower()

    # GAFF parameters from the downloadable package
    # prm_file = PACKAGE_DIR + "/gaff2.dat"

    amber = GeneralAmber(xyz_file=xyz_file)
    amber.compute_topology(periodic=periodic, simple=simple)
    if par_file is None:
        par_file = xyz_file.replace(".xyz", "(computed_parameters).par")
    amber.write_par(xyz_file.replace(".xyz", "(computed_parameters).par"))
    amber.atoms.write_xyz(xyz_file.replace(".xyz", ".types.xyz"))

    input("Press anything when the par file is ready.")

    lmp = LmpDat()
    lmp.atoms = amber.atoms
    lmp.get_params(par_file)
    lmp.cell_to_lo_hi()
    lmp.write_lmpdat(lmpdat_file)


def write_lmpdat_compass(xyz_file: str, lmpdat_file: str, par_file=None,
                         periodic="", simple=False):

    """
    UNDER DEVELOPMENT!
    """

    # gambiarra to fix my input xyz file
    pdms = True

    periodic = periodic.lower()

    compass = Compass(xyz_file=xyz_file)

    if pdms:
        for atom in compass.atoms:
            if atom.type == "S":
                atom.type = "Si"

    compass.compute_topology(periodic=periodic, simple=simple)
    #compass.order_in_bonds()
    if par_file is None:
        par_file = xyz_file.replace(".xyz", ".par")
    compass.write_par(xyz_file.replace(".xyz", ".par"))
    compass.atoms.write_xyz(xyz_file.replace(".xyz", ".types.xyz"),
                            with_classification=True)

    input("Press anything when the par file is ready.")

    lmp = LmpDat()
    lmp.atoms = compass.atoms
    compass.order_in_bonds()
    lmp.get_params(par_file)

    # COMPASS-specific adjustments; all are necessary!
    compass.order_in_bonds()
    compass.increment_charges()
    compass.order_in_angles()
    compass.order_in_dihedrals()

    lmp.cell_to_lo_hi()
    lmp.write_lmpdat(lmpdat_file)


def finish_lmpdat(xyz_types_file: str, lmpdat_file: str, par_file: str,
                  periodic="", simple=False):
    """
    Writes a LAMMPS Data file with atom_style='full' from classified xyz file,
    which can be generated with the write_lmpdat functions.

    This function (finish_lmpdat) is supposed to be used when a write_lmpdat
    stops before finishing (if the par file is not ready when it should be),
    or if there is an xyz file with the wanted atom types, as well as a par
    file with the corresponding parameters (e.g. for a force field other than
    the ones implemented in this package).

    Parameters
    ----------
    xyz_types_file : str
        Path to the input xyz file previously generated by the
        write_lmpdat function (atoms already classified).
    lmpdat_file : str
        Path to the output LAMMPS Data file.
    par_file : str
        Path to the input parameters file (.par) previously generated.
    periodic : str, optional
        Axes in which the system is periodic (e.g. 'xyz', 'xy', 'z';
        use an empty string '' for non-periodic). Standard is ''.
    simple : bool, optional
        If a simple, less optimised algorithm is wanted for computing
        bonds, angles and dihedrals. (For small systems. See Notes.)
        Standard is False.

    Notes
    -----
    If a system is periodic in all directions, but one of its dimensions is
    such that it's only divided into 1 or 2 regions, simple=True should be
    used to avoid bugs.

    Examples
    --------
    See directory 'examples'.

    """

    lmp = LmpDat()
    lmp.get_xyz(xyz_types_file)
    lmp.atoms.compute_topology(periodic=periodic, simple=simple)
    lmp.get_params(par_file)
    lmp.write_lmpdat(lmpdat_file)
