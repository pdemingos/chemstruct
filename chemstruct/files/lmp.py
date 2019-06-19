# Copyright 2019 Pedro G. Demingos

"""
Defines classes for dealing with LAMMPS input and output files.

The LmpDat class is for LAMMPS Data files, which contain
atomic coordinates and possibly force field parameters for
molecular simulation.

LAMMPS website: https://lammps.sandia.gov
LAMMPS documentation: https://lammps.sandia.gov/doc/Manual.html#

"""

import os
import numpy as np

from atoms import (Atom, Atoms, AtomType, BondType, AngleType, DihedralType,
                   ImproperType)
from tools import linear_counter
from files.main import File, find_between, clear_end
from files.csv import Csv
from files.xyz import Xyz


class LmpDat(File):
    """Class for LAMMPS Data files.
    Contains topological information for input.
    May contain force field parameters."""

    def __init__(self, absolute_path=None):
        super().__init__(absolute_path)
        self.tags.add("lmpdat")
        self._multiple_dihedrals = False

        self.atoms = Atoms([])
        self._xyz_lo_hi = None  # [xlo, xhi, ylo, yhi, zlo, zhi]

        if self.absolute_path is not None:
            self.read_lmpdat()

    def read_lmpdat(self):
        pass

    def get_xyz(self, xyz_file: str):
        """
        Gets info from xyz file, and lattice if there is any.

        Parameters
        ----------
        xyz_file : str
            Path to xyz file.

        """

        # makes a dummy Xyz object
        _xyz = Xyz(xyz_file)

        # uses the dummy Xyz object to get info from xyz file
        if len(self.atoms) > 0:
            print("WARNING: erasing current {} atoms".format(len(self.atoms)))
        self.atoms = _xyz.atoms
        print("Reading {} atoms from xyz file".format(len(self.atoms)))

        # deals with the cell
        self.cell_to_lo_hi()

    def cell_to_lo_hi(self):
        xlo = min(atom.position[0] for atom in self.atoms)
        ylo = min(atom.position[1] for atom in self.atoms)
        zlo = min(atom.position[2] for atom in self.atoms)
        extra_spacing = 5  # angstroms
        xhi = max(atom.position[0] for atom in self.atoms) + extra_spacing
        yhi = max(atom.position[1] for atom in self.atoms) + extra_spacing
        zhi = max(atom.position[2] for atom in self.atoms) + extra_spacing
        try:  # if self.atoms.cell is not None
            xhi = xlo + self.atoms.cell[0][0]
            yhi = ylo + self.atoms.cell[1][1]
            zhi = zlo + self.atoms.cell[2][2]
        except TypeError:  # if self.atoms.cell is None
            print("WARNING: no Lattice found in xyz file; "
                  "if the system is periodic, this is BAD")
        finally:
            self._xyz_lo_hi = [xlo, xhi, ylo, yhi, zlo, zhi]

    def get_params(self, params_file: str):
        """
        Gets info from parameters file (force field parameters).

        Parameters
        ----------
        params_file : str
            Path to parameters file.

        Raises
        ------
        TypeError
            If any parameter in the file isn't a number.

        Notes
        -----
        The expected file format isn't native to LAMMPS.
        Instead, it can be written by e.g. a CharmmGeneral object.
        The file should look like this:

        "
        Atom Types

        C1	-0.18	0.06	3.59923	0.01	3.38542
        C2	0.0	    0.06	3.59923	0.01	3.38542

        Bond Types

        C1:C1	195.0	1.53
        C1:C2	195.0	1.53

        Angle Types

        C1:C2:C1	 58.0	109.5	11.16	2.561
        C1:C1:C2     35.0	111.4	22.53	2.179

        Dihedral Types

        C1:C2:C1:C1	0.14	3	0	1.0

        Improper Types

        # none
        "

        """

        if not self.atoms.bonds:
            print("WARNING: topology not computed yet")

        # makes a dummy generic File for parameters
        _params_file = File(params_file)

        atom_types_index = _params_file._find("Atom Types")
        bond_types_index = _params_file._find("Bond Types")
        angle_types_index = _params_file._find("Angle Types")
        dihedral_types_index = _params_file._find("Dihedral Types")
        improper_types_index = _params_file._find("Improper Types")

        # each of the titles above must appear exactly once
        # note that the area below the title may be left empty

        for lis in [atom_types_index, bond_types_index,
                    angle_types_index, dihedral_types_index,
                    improper_types_index]:
            if (len(lis) > 1) or (len(lis) == 0):
                print("WARNING: bad params file")

        atom_types_index = atom_types_index[0]
        bond_types_index = bond_types_index[0]
        angle_types_index = angle_types_index[0]
        dihedral_types_index = dihedral_types_index[0]
        improper_types_index = improper_types_index[0]

        # checks if there are multiple dihedrals
        self._multiple_dihedrals = False
        types = []
        for line in _params_file.content[dihedral_types_index + 1:improper_types_index]:
            try:
                if line[0].isalpha():
                    typ, *etc = tuple(line.split())
                    if typ in types:
                        self._multiple_dihedrals = True
                    types.append(typ)
            except IndexError:
                continue

        # Atom Types

        index = 1
        for line in _params_file.content[atom_types_index + 1:bond_types_index]:

            try:
                if line[0].isalpha():
                    typ, *charge_and_params = tuple(line.split())
                    if typ not in AtomType.instances_dict.keys():
                        AtomType(typ)  # instantiates it
                    charge, *params = charge_and_params
                    params_nums = []

                    # just takes strings into the right numeric types
                    # if the str "1" becomes the float 1.0, LAMMPS may complain
                    for param in params:
                        try:
                            if "." in param:
                                param_num = float(param)
                            else:
                                param_num = int(param)
                        except TypeError:
                            raise TypeError("all Atom Type parameters must be numbers, "
                                            "at least one is not")
                        except ValueError:  # if param=None
                            param_num = None
                        params_nums.append(param_num)

                    AtomType.instances_dict[typ].index = index
                    AtomType.instances_dict[typ].charge = charge
                    AtomType.instances_dict[typ].params = params_nums  # list
                    index += 1
            except IndexError:
                continue

        # Bond Types

        index = 1
        for line in _params_file.content[bond_types_index + 1:angle_types_index]:

            try:
                if line[0].isalpha():
                    typ, *params = tuple(line.split())
                    if typ not in BondType.instances_dict.keys():
                        BondType(typ)  # instantiates it
                    params_nums = []

                    # just takes strings into the right numeric types
                    # if the str "1" becomes the float 1.0, LAMMPS may complain
                    for param in params:
                        try:
                            if "." in param:
                                param_num = float(param)
                            else:
                                param_num = int(param)
                        except TypeError:
                            raise TypeError("all Bond Type parameters must be numbers, "
                                            "at least one is not")
                        except ValueError:  # if param is None
                            param_num = None
                        params_nums.append(param_num)

                    BondType.instances_dict[typ].index = index
                    BondType.instances_dict[typ].params = params_nums  # list
                    index += 1
            except IndexError:
                continue

        # Angle Types

        index = 1
        for line in _params_file.content[angle_types_index + 1:dihedral_types_index]:

            try:
                if line[0].isalpha():
                    typ, *params = tuple(line.split())
                    if typ not in AngleType.instances_dict.keys():
                        AngleType(typ)  # instantiates it
                    params_nums = []

                    # just takes strings into the right numeric types
                    # if the str "1" becomes the float 1.0, LAMMPS may complain
                    for param in params:
                        try:
                            if "." in param:
                                param_num = float(param)
                            else:
                                param_num = int(param)
                        except TypeError:
                            raise TypeError("all Angle Type parameters must be numbers, "
                                            "at least one is not")
                        except ValueError:  # if param=None
                            param_num = None
                        params_nums.append(param_num)

                    AngleType.instances_dict[typ].index = index
                    AngleType.instances_dict[typ].params = params_nums  # list
                    index += 1
            except IndexError:
                continue

        # Dihedral Types

        index = 1
        for line in _params_file.content[dihedral_types_index + 1:improper_types_index]:

            try:
                if line[0].isalpha():
                    typ, *params = tuple(line.split())
                    if typ not in DihedralType.instances_dict.keys():
                        DihedralType(typ)  # instantiates it
                    params_nums = []

                    # just takes strings into the right numeric types
                    # if the str "1" becomes the float 1.0, LAMMPS may complain
                    for param in params:
                        try:
                            if "." in param:
                                param_num = float(param)
                            else:
                                param_num = int(param)
                        except TypeError:
                            raise TypeError("all Dihedral Type parameters must be numbers, "
                                            "at least one is not")
                        except ValueError:  # if param=None
                            param_num = None
                        params_nums.append(param_num)

                    if not self._multiple_dihedrals:
                        DihedralType.instances_dict[typ].index = index
                        DihedralType.instances_dict[typ].params = params_nums  # list
                    else:  # if self._multiple_dihedrals
                        try:
                            DihedralType.instances_dict[typ].index.append(index)
                            DihedralType.instances_dict[typ].params.append(params_nums)
                        except AttributeError:  # first params of each dihedral
                            DihedralType.instances_dict[typ].index = []
                            DihedralType.instances_dict[typ].params = []
                            DihedralType.instances_dict[typ].index.append(index)
                            DihedralType.instances_dict[typ].params.append(params_nums)
                    index += 1
            except IndexError:
                continue

        # Improper Types

        index = 1
        for line in _params_file.content[improper_types_index + 1:]:

            try:
                if line[0].isalpha():
                    typ, *params = tuple(line.split())
                    if typ not in ImproperType.instances_dict.keys():
                        ImproperType(typ)
                    params_nums = []

                    # just takes strings into the right numeric types
                    # if the str "1" becomes the float 1.0, LAMMPS may complain
                    for param in params:
                        try:
                            if "." in param:
                                param_num = float(param)
                            else:
                                param_num = int(param)
                        except TypeError:
                            raise TypeError("all Improper Type parameters must be numbers, "
                                            "at least one is not")
                        except ValueError:  # if param=None
                            param_num = None
                        params_nums.append(param_num)

                    ImproperType.instances_dict[typ].index = index
                    ImproperType.instances_dict[typ].params = params_nums  # list
                    index += 1
            except IndexError:
                continue

        # moves charge from AtomType to each Atom
        self._set_charges()

        # re-indexes types leaving extra types (in .par but not in .xyz) behind
        self.atoms.re_index_types()

        # GAMBIARRA
        # but should not be a problem if .par is right
        self.delete_impropers_without_parameters()

    def delete_impropers_without_parameters(self):

        impropers = self.atoms.impropers[:]
        i = 1
        for improper in impropers:
            if improper.type.params is None:
                self.atoms.impropers.remove(improper)
            else:
                improper.index = i
                i += 1

        improper_types = self.atoms.improper_types[:]
        i = 1
        for improper_type in improper_types:
            if improper_type.params is None:
                self.atoms.improper_types.remove(improper_type)
            else:
                improper_type.index = i
                i += 1

    def _set_charges(self):
        for atom in self.atoms:
            atom.charge = atom.type.charge

    def multiply_dihedral_coeffs(self):
        """
        Turns every dihedral in n dihedrals, where n is the number of
        terms in its series (e.g. a Fourier Series).

        Returns
        -------
        n_dihedrals : int
            Total number of dihedrals in self.atoms, computed
            after the multiplication.

        Notes
        -----
        The CHARMM parameters are given in this format.
        LAMMPS also has a nice sinusoidal form for dihedrals,
        making this format well suited.

        """

        if not self._multiple_dihedrals:
            raise TypeError("not self._multiple_dihedrals")
            # this means there's an internal problem

        for dihedral in self.atoms.dihedrals:
            try:
                dihedral.index = dihedral.type.index[:]
                # meaningless numbers, only the list's len matters
            except TypeError:
                raise TypeError("bad dihedral type {}, check .par file".format(
                    str(dihedral.type)))

        index = 1
        for dihedral in self.atoms.dihedrals:
            for i in range(len(dihedral.index)):
                dihedral.index[i] = index
                index += 1

        n_dihedrals = index - 1
        return n_dihedrals

    def write_lmpdat(self, filename, atom_style="full", comments=False):
        """
        Writes LAMMPS Data file with info present in this object
        i.e. topology and force field parameters previously read.

        Parameters
        ----------
        filename : str
            Path for output LAMMPS Data file.
            CARE: Will erase any existing file with this name.
        atom_style : {'full', ...}, optional
            Wanted atom_style output format.
        comments : bool, optional
            If comments are wanted in the output file,
            containing redundant info about atom types.

        Notes
        -----
        See LAMMPS documentation on atom_style:
        https://lammps.sandia.gov/doc/atom_style.html

        Examples
        --------
        Example for combining xyz file and par file to make lmpdat file:
        >> xyz = "/home/lasim/mydir/NAME.xyz"
        >> par = "/home/lasim/mydir/NAME.par"
        >> from files.lmp import LmpDat
        >> lmp = LmpDat()
        >> lmp.get_xyz(xyz)
        >> lmp.atoms.compute_topology()  # finds bonds, angles, etc
        >> lmp.get_params(par)  # AFTER computing topology
        >> lmp.write_lmpdat("/home/lasim/mydir/NAME.lmp")

        """

        # sets flags according to atom_style
        if atom_style == "full":
            atoms, bonds, angles, dihedrals, impropers = True, True, True, True, True
            molecule, charge = True, True
            parameters = True
        elif atom_style == "atomic":
            atoms, bonds, angles, dihedrals, impropers = True, False, False, False, False
            molecule, charge = False, False
            parameters = False
        elif atom_style == "charge":
            atoms, bonds, angles, dihedrals, impropers = True, False, False, False, False
            molecule, charge = False, True
            parameters = False
        else:  # to be expanded as needed
            raise TypeError("bad atom_style")

        number_of_dihedrals = len(self.atoms.dihedrals)
        if self._multiple_dihedrals:
            number_of_dihedrals = self.multiply_dihedral_coeffs()

        with open(filename, "w") as F:

            # HEADERS

            F.write("LAMMPS data file\n\n")

            if atoms:
                F.write(str(len(self.atoms.atoms)) + " atoms \n")
            if bonds:
                F.write(str(len(self.atoms.bonds)) + " bonds \n")
            if angles:
                F.write(str(len(self.atoms.angles)) + " angles \n")
            if dihedrals:
                F.write(str(number_of_dihedrals) + " dihedrals \n")
            if impropers:
                F.write(str(len(self.atoms.impropers)) + " impropers \n\n")

            if atoms:
                F.write(str(len(self.atoms.atom_types)) + " atom types \n")
            if bonds:
                F.write(str(len(self.atoms.bond_types)) + " bond types \n")
            if angles:
                F.write(str(len(self.atoms.angle_types)) + " angle types \n")
            if dihedrals and not self._multiple_dihedrals:
                F.write(str(len(self.atoms.dihedral_types)) + " dihedral types \n")
            elif dihedrals and self._multiple_dihedrals:
                F.write(str(sum(len(typ.index) for typ in self.atoms.dihedral_types)) + " dihedral types \n")
            if impropers:
                F.write(str(len(self.atoms.improper_types)) + " improper types \n\n")

            F.write(str(self._xyz_lo_hi[0]) + " " + str(self._xyz_lo_hi[1]) + " " + "xlo xhi \n")
            F.write(str(self._xyz_lo_hi[2]) + " " + str(self._xyz_lo_hi[3]) + " " + "ylo yhi \n")
            F.write(str(self._xyz_lo_hi[4]) + " " + str(self._xyz_lo_hi[5]) + " " + "zlo zhi \n\n")

            F.write("Masses \n\n")
            for typ in self.atoms.atom_types:
                F.write(str(typ.index) + " " + str(typ.mass) + "  # " + str(typ) + "\n")

            # STRUCTURES

            if atoms and self.atoms.atoms:
                molecule_flag = False
                F.write("\nAtoms \n\n")
                if molecule and charge:
                    for atom in self.atoms.atoms:
                        if atom.molecule.index is None:
                            atom.molecule.index = 1
                            molecule_flag = True
                        F.write(str(atom.index) + " " + str(atom.molecule.index) + " " + str(atom.type.index) +
                                " " + str(atom.charge) + " " + " ".join(str(pos) for pos in atom.position)
                                + "  # " + str(atom.type) + "\n")
                elif charge:
                    for atom in self.atoms.atoms:
                        if atom.charge is None:
                            atom.charge = 0.0
                        F.write(str(atom.index) + " " + str(atom.type.index) +
                                " " + str(atom.charge) + " " + " ".join(str(pos) for pos in atom.position)
                                + "  # " + str(atom.type) + "\n")
                else:
                    for atom in self.atoms.atoms:
                        F.write(str(atom.index) + " " + str(atom.type.index) +
                                " " + " ".join(str(pos) for pos in atom.position)
                                + "  # " + str(atom.type) + "\n")

                if molecule_flag:
                    print("WARNING: single molecule assumed!")

            if bonds and self.atoms.bonds:
                F.write("\nBonds \n\n")
                for bond in self.atoms.bonds:
                    if comments:
                        output = (str(bond.index) + " " + str(bond.type.index) + " " +
                                  " ".join(str(atom.index) for atom in bond.atoms)
                                  + "  # " + str(bond.type) + "\n")
                    else:
                        output = (str(bond.index) + " " + str(bond.type.index) + " " +
                                  " ".join(str(atom.index) for atom in bond.atoms) + "\n")
                    F.write(output)

            if angles and self.atoms.angles:
                F.write("\nAngles \n\n")
                for angle in self.atoms.angles:
                    if comments:
                        output = (str(angle.index) + " " + str(angle.type.index) + " " +
                                  " ".join(str(atom.index) for atom in angle.atoms)
                                  + "  # " + str(angle.type) + "\n")
                    else:
                        output = (str(angle.index) + " " + str(angle.type.index) + " " +
                                  " ".join(str(atom.index) for atom in angle.atoms) + "\n")
                    F.write(output)

            if dihedrals and self.atoms.dihedrals:
                F.write("\nDihedrals \n\n")
                for dihedral in self.atoms.dihedrals:
                    if not self._multiple_dihedrals:
                        if comments:
                            output = (str(dihedral.index) + " " + str(dihedral.type.index) + " " +
                                      " ".join(str(atom.index) for atom in dihedral.atoms)
                                      + "  # " + str(dihedral.type) + "\n")
                        else:
                            output = (str(dihedral.index) + " " + str(dihedral.type.index) + " " +
                                      " ".join(str(atom.index) for atom in dihedral.atoms) + "\n")
                        F.write(output)
                    else:  # if self._multiple_dihedrals
                        for i in range(len(dihedral.index)):
                            if comments:
                                output = (str(dihedral.index[i]) + " " + str(dihedral.type.index[i]) + " " +
                                          " ".join(str(atom.index) for atom in dihedral.atoms)
                                          + "  # " + str(dihedral.type) + "\n")
                            else:
                                output = (str(dihedral.index[i]) + " " + str(dihedral.type.index[i]) + " " +
                                          " ".join(str(atom.index) for atom in dihedral.atoms) + "\n")
                            F.write(output)

            if impropers and self.atoms.impropers:
                F.write("\nImpropers \n\n")
                for improper in self.atoms.impropers:
                    if comments:
                        output = (str(improper.index) + " " + str(improper.type.index) + " " +
                                  " ".join(str(atom.index) for atom in improper.atoms)
                                  + "  # " + str(improper.type) + "\n")
                    else:
                        output = (str(improper.index) + " " + str(improper.type.index) + " " +
                                  " ".join(str(atom.index) for atom in improper.atoms) + "\n")
                    F.write(output)

            # TYPES

            if parameters and atoms and self.atoms.atoms:
                F.write("\nPair Coeffs \n\n")
                for atom_type in self.atoms.atom_types:
                    F.write(str(atom_type.index) + " " +
                            " ".join(str(param) for param in atom_type.params) +
                            "  # " + str(atom_type) + "\n")

            if bonds and self.atoms.bonds:
                F.write("\nBond Coeffs \n\n")
                for bond_type in self.atoms.bond_types:
                    F.write(str(bond_type.index) + " " +
                            " ".join(str(param) for param in bond_type.params) +
                            "  # " + str(bond_type) + "\n")

            if angles and self.atoms.angles:
                F.write("\nAngle Coeffs \n\n")
                for angle_type in self.atoms.angle_types:
                    F.write(str(angle_type.index) + " " +
                            " ".join(str(param) for param in angle_type.params) +
                            "  # " + str(angle_type) + "\n")

            if dihedrals and self.atoms.dihedrals:
                F.write("\nDihedral Coeffs \n\n")
                for dihedral_type in self.atoms.dihedral_types:
                    if not self._multiple_dihedrals:
                        F.write(str(dihedral_type.index) + " " +
                                " ".join(str(param) for param in dihedral_type.params) +
                                "  # " + str(dihedral_type) + "\n")
                    else:  # if self._multiple_dihedrals
                        for i in range(len(dihedral_type.index)):
                            F.write(str(dihedral_type.index[i]) + " " +
                                    " ".join(str(param) for param in dihedral_type.params[i]) +
                                    "  # " + str(dihedral_type) + "\n")

            if impropers and self.atoms.impropers:
                F.write("\nImproper Coeffs \n\n")
                for improper_type in self.atoms.improper_types:
                    try:
                        F.write(str(improper_type.index) + " " +
                                " ".join(str(param) for param in improper_type.params) +
                                "  # " + str(improper_type) + "\n")
                    except TypeError:
                        print("WARNING: improper {} without parameters!".format(
                            str(improper_type)))


class LmpOut(File):
    """Class for LAMMPS output files."""

    def __init__(self, absolute_path):
        super().__init__(absolute_path=absolute_path)

        self.tables = []  # list of Csv objects; each Csv is a table

        if self.absolute_path is not None:
            self.read()

    def read(self):
        start_indices = self._find("Per MPI rank memory allocation")
        end_indices = self._find("Loop time")

        # can add more info to compute

        for start_and_end_indices in zip(start_indices, end_indices):
            start_i, end_i = start_and_end_indices
            csv = Csv()

            headers = tuple(self.content[start_i + 1].split())
            columns = [[] for _ in range(len(headers))]
            for line in self.content[start_i + 2:end_i]:
                terms = tuple(line.split())
                for i in range(len(headers)):
                    columns[i].append(float(terms[i]))

            for i in range(len(headers)):
                csv.add_header(headers[i])
                csv.add_array(columns[i])
            self.tables.append(csv)

    def write_tables(self, path=None):
        """
        Writes a csv file for each table computed from the LAMMPS output.
        Files are numbered in the order they were read from the output.

        Parameters
        ----------
        path : str, optional
            Path to the csv file. If not given, it's written in the same
            directory as the LAMMPS output file.

        """
        if path is None:
            path = self.absolute_path.replace(".out", "")
        else:
            path = path.replace(".csv", "")

        for (i, csv) in enumerate(self.tables):
            this_path = path + "({}).csv".format(i + 1)
            csv.write_csv(this_path)

    def write_single_table(self, path=None):
        """
        Writes a single csv file for all tables computed from the output.

        Parameters
        ----------
        path : str, optional
            Path to the csv file. If not given, it's written in the same
            directory as the LAMMPS output file.

        """
        if path is None:
            path = self.absolute_path.replace(".out", ".csv")

        single_table = Csv()
        for header in self.tables[0].headers:
            single_table.add_array([])
            single_table.add_header(header)
        for csv in self.tables:
            for (i, column) in enumerate(csv.data):
                single_table.data[i] += column
        single_table.max_len = len(single_table.data[0])
        single_table.write_csv(path)


class Cfg(File):
    """Class for LAMMPS cfg files."""
    # usual line: "mass type x y z *etc" where len(x,y,z,*etc) == entry_count

    def __init__(self, absolute_path=None,
                 ignore_types=None):  # ignore_types is for internal use
        super().__init__(absolute_path=absolute_path)

        self.number_of_particles = None
        self.has_velocity = True
        self.atoms = Atoms([])
        self.entry_count = None
        self.auxiliary = []

        if ignore_types is None:
            self.ignore_types = []
        else:
            self.ignore_types = ignore_types

        if self.absolute_path is not None:
            self.read()

    def read(self):

        index_start = self._find("Number of particles")[0]
        index_auxiliary = None
        index_atoms = None
        self.number_of_particles = int(self.content[index_start].split()[-1])

        h0_11, h0_12, h0_13, h0_21, h0_22, h0_23, h0_31, h0_32, h0_33 = 0, 0, 0, 0, 0, 0, 0, 0, 0

        if len(self._find(".NO_VELOCITY.")) == 1:
            self.has_velocity = False

        if len(self._find("A = 1 Angstrom")) == 0:
            print("WARNING: distance units in cfg file may not be in Angstroms!")

        for (index, line) in enumerate(self.content[index_start:], index_start):
            if line.startswith("H0"):
                if "H0(1,1)" in line:
                    h0_11 = float(line.split()[-2])
                elif "H0(1,2)" in line:
                    h0_12 = float(line.split()[-2])
                elif "H0(1,3)" in line:
                    h0_13 = float(line.split()[-2])
                elif "H0(2,1)" in line:
                    h0_21 = float(line.split()[-2])
                elif "H0(2,2)" in line:
                    h0_22 = float(line.split()[-2])
                elif "H0(2,3)" in line:
                    h0_23 = float(line.split()[-2])
                elif "H0(3,1)" in line:
                    h0_31 = float(line.split()[-2])
                elif "H0(3,2)" in line:
                    h0_32 = float(line.split()[-2])
                elif "H0(3,3)" in line:
                    h0_33 = float(line.split()[-2])
            elif "entry_count" in line:
                self.entry_count = int(line.split()[-1])
                n_standard_args = 6 if self.has_velocity else 3
                for _ in range(self.entry_count - n_standard_args):
                    self.auxiliary.append(None)
                index_auxiliary = index + 1
                break
            else:
                continue

        self.atoms.cell = [[h0_11, h0_12, h0_13],
                           [h0_21, h0_22, h0_23],
                           [h0_31, h0_32, h0_33]]

        for (index, line) in enumerate(self.content[index_auxiliary:], index_auxiliary):
            if line.startswith("auxiliary"):
                i = int(find_between(line, "[", "]"))
                _, __, auxiliary = tuple(line.split())
                self.auxiliary[i] = auxiliary
            else:
                index_atoms = index
                break

        for i in range(index_atoms, len(self.content), 3):
            typ = clear_end(self.content[i + 1], [" ", "\n", "\t"])
            if typ in self.ignore_types:
                continue
            mass = float(self.content[i + 0])
            etc = dict()
            if self.has_velocity:
                raise TypeError("cfg with velocities! not implemented yet!")
            else:
                x, y, z, *args = tuple(self.content[i + 2].split())
                atom = Atom(atom_type=typ, position=[float(x), float(y), float(z)])
                atom.position = np.matmul(self.atoms.cell, atom.position)
                atom.type.mass = mass
                for (index, arg) in enumerate(args):
                    etc[self.auxiliary[index]] = arg  # string
                atom.etc = etc
                self.atoms.add_atom(atom)

        if self.number_of_particles != len(self.atoms) and not self.ignore_types:
            print("WARNING: file says {} atoms, but only {} were found".format(
                self.number_of_particles, len(self.atoms)))

    def write_xyz(self, path, real_types=False):
        xyz = Xyz()
        xyz.atoms = self.atoms
        xyz.write_xyz(path, real_types=real_types)
