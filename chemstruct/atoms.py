# Copyright 2019 Pedro G. Demingos

"""
Defines Atom and Atoms classes, as well as
topological structures and topological types.

"""

import numpy as np
import math
from time import time

from constants import BOND_LENGTHS, ATOMIC_MASSES, axis_to_dim
from tools import break_regions


class Atom:
    """Class for an Atom."""

    index = 1

    def __init__(self, **kwargs):

        self.index = Atom.index
        Atom.index += 1
        # every Atom object is instantiated with a unique index
        # but within an Atoms object it may have its index changed

        self._type = None
        self.type = kwargs.get("atom_type", "Atom")
        self._position = None
        self.position = kwargs.get("position", [0, 0, 0])
        self.if_pos = kwargs.get("if_pos", [1, 1, 1])  # QE
        self.charge = kwargs.get("charge", None)
        self.etc = dict()

        # molecular info
        self.neighbors = set()  # set of atoms this atom is bonded to
        self.molecule = None
        self.hybridization = None
        self.cycles = []
        self.topological_tags = set()
        self.classification = None

    def __str__(self, real_type=False):
        return self._type.__str__(real_type=real_type) + \
               "\t" + "\t".join([str(p) for p in self._position])

    def __repr__(self):
        return self._type.__str__()

    def copy(self):
        """
        Copies the Atom object, instantiating a new object with the same
        type, position, if_pos and charge. Doesn't change the original.

        Returns
        -------
        copy : Atom
            The atom's copy.

        """
        copy = Atom(atom_type=self.type,
                    position=self.position,
                    if_pos=self.if_pos,
                    charge=self.charge)  # add more if necessary
        # carries no neighbors, hybridization, etc
        return copy

    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, atom_type):
        if isinstance(atom_type, AtomType):
            self._type = atom_type
        elif isinstance(atom_type, str):
            try:  # if the AtomType already exists
                self._type = AtomType.instances_dict[atom_type]
            except KeyError:  # instantiates new AtomType
                self._type = AtomType(atom_type)
        else:
            raise TypeError("atom type must be a string or an AtomType object,"
                            " got type {}".format(type(atom_type)))

    @property
    def position(self):
        return self._position

    @position.setter
    def position(self, position):
        try:
            self._position = np.array(position, dtype=float)
        except TypeError:
            raise TypeError("bad position type, got {}".format(type(position)))
        except ValueError:
            raise ValueError("items in position must be numbers")
        self._position = np.round(self._position, 6)

    def get_neighbors(self, neighbor_type_real: str):
        """Returns a list of all the atom's neighbors with given real type.
        For example, the real type of some carbon 'C3' is 'C'."""
        others = []
        for other in self.neighbors:
            if other.type.real == neighbor_type_real:
                others.append(other)
        return others

    def fix(self, axis: str):
        """Sets the atom's if_pos to 0 in the given axis.
        Used for QUANTUM ESPRESSO interface only."""
        dim = axis_to_dim(axis)
        self.if_pos[dim] = 0

    def translate(self, vector):
        """Translates (moves) the atom's position by the given vector."""
        self.position = self.position + np.array(vector)

    def rotate(self, angle: float, axis: str):
        """
        Rotates atom by the given angle around the given axis.

        Parameters
        ----------
        angle : float
            Angle (in degrees) of rotation.
        axis : str
            Axis to rotate around: 'x', 'y' or 'z'.

        """
        angle *= np.pi / 180.0  # deg to rad
        sin = np.sin(angle)
        cos = np.cos(angle)
        axis = axis.lower()
        if axis == "x":
            rotation_matrix = [[1, 0, 0], [0, cos, -sin], [0, sin, cos]]
        elif axis == "y":
            rotation_matrix = [[cos, 0, sin], [0, 1, 0], [-sin, 0, cos]]
        elif axis == "z":
            rotation_matrix = [[cos, -sin, 0], [sin, cos, 0], [0, 0, 1]]
        else:
            raise ValueError("axis should be 'x', 'y' or 'z', got {}".format(
                axis))
        rotation_matrix = np.array(rotation_matrix)
        self.position = np.matmul(rotation_matrix, self.position)


class Atoms:
    """Class for a system of atoms.
    May contain topological info (bonds, angles etc)."""

    def __init__(self, atoms: list):

        # topological structures
        self.atoms = atoms
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.impropers = []
        self.cycles = [None, None, None, []]  # 0, 1, 2, 3
        self.molecules = []

        # topological types
        self.atom_types = []
        self.bond_types = []
        self.angle_types = []
        self.dihedral_types = []
        self.improper_types = []

        self.tags = set()
        self._cell = None

    def __len__(self):
        return len(self.atoms)

    def __iter__(self):
        return iter(self.atoms)

    def __getitem__(self, index_or_slice):
        # taking a slice of atoms doesn't carry topological structures
        return self.atoms[index_or_slice]  # Atom or list

    def __str__(self):
        print_atoms = ""
        for atom in self.atoms:
            print_atoms += atom.__str__() + "\n"
        return print_atoms

    def __repr__(self):
        return "Atoms( n_atoms={} )".format(len(self.atoms))

    def __add__(self, other):
        # does NOT keep the cycles!
        if isinstance(other, Atoms):
            new = Atoms([atom for atom in self.atoms + other.atoms])

            # topological structures
            new.bonds = self.bonds + other.bonds
            new.angles = self.angles + other.angles
            new.dihedrals = self.dihedrals + other.dihedrals
            new.impropers = self.impropers + other.impropers
            new.molecules = self.molecules + other.molecules

            # topological types
            new.atom_types = self.atom_types
            new.atom_types += [typ for typ in other.atom_types
                               if typ not in self.atom_types]
            new.bond_types = self.bond_types
            new.atom_types += [typ for typ in other.bond_types
                               if typ not in self.bond_types]
            new.angle_types = self.angle_types
            new.atom_types += [typ for typ in other.angle_types
                               if typ not in self.angle_types]
            new.dihedral_types = self.dihedral_types
            new.atom_types += [typ for typ in other.dihedral_types
                               if typ not in self.dihedral_types]
            new.improper_types = self.improper_types
            new.atom_types += [typ for typ in other.improper_types
                               if typ not in self.improper_types]
            new.re_index_types()
            return new

        else:
            return NotImplemented

    @property
    def cell(self):
        return self._cell

    @cell.setter
    def cell(self, cell):
        if (not isinstance(cell, list)) and (not isinstance(cell, np.ndarray)):
            raise TypeError("Bad cell type: expected list or numpy array,"
                            " got {}".format(type(cell)))
        try:
            dummy_cell = np.array(cell, dtype=float)
        except ValueError:
            raise ValueError("Bad cell shape or content: expected shape=(3,3),"
                             " filled with only numbers")
        if dummy_cell.shape != (3, 3):
            raise ValueError("Bad cell shape: expected shape=(3,3),"
                             " got {}".format(dummy_cell.shape))
        self._cell = dummy_cell

    def copy(self, topology=False):
        """Returns new Atoms object with the same topological structures."""
        copy = Atoms([])
        atom_map = dict()  # for copying topology
        for atom in self.atoms:
            atom_copy = atom.copy()
            copy.add_atom(atom_copy)
            atom_map[atom.index] = atom_copy
        if self.cell is not None:
            copy.cell = self.cell
        if topology:
            for bond in self.bonds:
                copy.add_bond([atom_map[bond[0].index],
                               atom_map[bond[1].index]])
            for angle in self.angles:
                copy.add_angle([atom_map[angle[0].index],
                                atom_map[angle[1].index],
                                atom_map[angle[2].index]])
            for dihedral in self.dihedrals:
                copy.add_dihedral([atom_map[dihedral[0].index],
                                   atom_map[dihedral[1].index],
                                   atom_map[dihedral[2].index],
                                   atom_map[dihedral[3].index]])
            for improper in self.impropers:
                copy.add_improper([atom_map[improper[0].index],
                                   atom_map[improper[1].index],
                                   atom_map[improper[2].index],
                                   atom_map[improper[3].index]])
            for (order, cycles_of_some_order) in enumerate(self.cycles):
                if cycles_of_some_order is None:
                    continue
                for cycle in cycles_of_some_order:
                    copy.add_cycle([atom_map[cycle[i].index]
                                    for i in range(order)])
            for molecule in self.molecules:
                copy.add_molecule([atom_map[atom.index]
                                   for atom in molecule.atoms])
        return copy

    def add_bond(self, bond: list):
        self.bonds.append(Bond(bond))
        bond[0].neighbors.add(bond[1])
        bond[1].neighbors.add(bond[0])

    def add_angle(self, angle: list):
        self.angles.append(Angle(angle))

    def add_dihedral(self, dihedral: list):
        dihedral = Dihedral(dihedral)
        if not dihedral.is_small_cycle:
            self.dihedrals.append(dihedral)

    def add_improper(self, improper: list):
        self.impropers.append(Improper(improper))

    def add_cycle(self, cycle: list):
        order = len(cycle)
        if len(set(cycle)) == order:  # else: repeated atom(s)
            while len(self.cycles) < order + 1:
                self.cycles.append([])
            if cycle not in self.cycles[order]:
                cycle_obj = Cycle(cycle)
                self.cycles[order].append(cycle_obj)
                # Atoms.cycles is a list of (Nones and) lists of cycle objects
                for atom in cycle:
                    atom.topological_tags.add("cycle {}".format(order))
                    atom.cycles.append(cycle_obj)
                    # Atom.cycles is a list of cycle objects

    def add_molecule(self, molecule: list):
        molecule = Molecule(molecule)
        self.molecules.append(molecule)
        for atom in molecule.atoms:
            atom.molecule = molecule

    def has_same_bonds(self, other_atoms):
        """Checks if two Atoms objects have the same bonds.
        The order of the atoms must be the same."""
        if len(self) != len(other_atoms):
            return False
        if len(self.bonds) != len(other_atoms.bonds):
            return False
        for (i, atom) in enumerate(self.atoms):
            other = other_atoms.atoms[i]
            # print("{}={}".format(i, atom.index))
            atom_neighbors = {n.index for n in atom.neighbors}
            other_neighbors = {n.index for n in other.neighbors}
            # print(atom_neighbors, other_neighbors)
            if atom_neighbors == other_neighbors:
                continue
            else:
                return False
        return True

    def compute_topology(self, periodic="", impropers=True, complete=False,
                         molecules=True, hold_pool_top_types=False,
                         bonds=True, simple=False):
        """
        Finds topological structures: bonds, angles and dihedrals.
        If asked, classifies atoms, finding hybridizations, molecules, etc.

        Parameters
        ----------
        periodic : str, optional
            Axes in which the system is periodic (e.g. 'xyz', 'xy', 'z';
            use an empty string '' for non-periodic). Standard is ''.
        impropers : bool, optional
            If improper dihedrals are wanted. Standard is True.
        complete : bool, optional
            If a complete topological analysis is wanted. This includes
            atom classifications. Standard is False.
        molecules : bool, optional
            If molecules are wanted. (Useful to turn off if there are
            macromolecules.) Standard is True.
        hold_pool_top_types : bool, optional
            Doesn't pool topological types. (For internal use, when the
            pool is called later.) Standard is False.
        bonds : bool, optional
            If the bonds are to be computed. (For internal use, when the
            bonds are already given.) Standard is True.
        simple : bool, optional
            If a simple, less optimised algorithm is wanted for computing
            bonds, angles and dihedrals. (For small systems. See Notes.)
            Standard is False.

        Notes
        -----
        This may take a while for large systems.

        Works with periodic triclinic systems.

        If the system is 1D, use periodic='z' or whatever makes sense.
        If the system is 2D, use periodic='xy' or whatever makes sense.
        Not doing this may cause problems.

        If a system is too small (in any axis), use simple=True.
        This will replace the optimised (divides-into-sub-regions)
        algorithm with a simpler, slower (checks-all-pairs) one.
        Not doing this may cause problems.

        """

        # bonds come first
        if bonds:
            # self.compute_bonds(periodic=periodic, triclinic=triclinic)
            # self.compute_bonds_smart(periodic=periodic)
            self.compute_smart(struct="bonds", periodic=periodic,
                               simple=simple)

        # then stuff that depend on bonds
        print("Computing other topological structures...")
        self.compute_cycles()
        # self.compute_angles()
        self.compute_smart(struct="angles", periodic=periodic,
                           simple=simple)
        # self.compute_dihedrals()
        self.compute_smart(struct="dihedrals", periodic=periodic,
                           simple=simple)
        self.compute_hybridizations()  # must come before impropers
        if molecules:
            self.compute_molecules()

        # then stuff that depend on the stuff above
        if impropers:
            self.compute_impropers()
        print("Other topological structures computed")

        # for specific atom types, e.g. for force field classification
        if complete:
            self.find_planar_cycles()
            # self.compute_bond_orders()
            self.compute_atom_classification()
            # ...

        # then the pools
        if not hold_pool_top_types:
            self.pool_topological_types()
            # does the indexing and molecular pools

    def move_atoms_into_cell(self):
        """Moves every atom into the cell, taken as a box between [0,0,0] and
        [lx,ly,lz]. Only makes sense for orthogonal cells."""
        count = 0
        for atom in self.atoms:
            was_moved = False
            for dim in [0, 1, 2]:
                if atom.position[dim] > self.cell[dim][dim]:
                    atom.position[dim] -= self.cell[dim][dim]
                    was_moved = True
                elif atom.position[dim] < 0:
                    atom.position[dim] += self.cell[dim][dim]
                    was_moved = True
            if was_moved:
                count += 1
        if count:
            print("Moved {} atoms into the box".format(count))

    def compute_smart(self, struct, periodic="", simple=False):
        """
        Computes bonds, angles or dihedrals.

        Parameters
        ----------
        struct : str
            What topological structure is wanted. Must be either
            'bonds', 'angles' or 'dihedrals'.
        periodic : str, optional
            Axes in which the system is periodic (e.g. 'xyz', 'xy', 'z';
            use an empty string '' for non-periodic). Standard is ''.
        simple : bool, optional
            If a simple, less optimised algorithm is wanted for computing
            bonds, angles and dihedrals. (For small systems. See Notes.)
            Standard is False.

        Notes
        -----
        Computing bonds requires checking pairs of atoms.
        Computing angles requires checking pairs of bonds.
        Computing dihedrals requires checking pairs of angles.
        This method uses an optimised algorithm that divides the system into
        regions, placing each object in one of them. Then, checks each region
        for internal pairs, as well as neighboring regions for inter-region
        pairs. This will cause problems for small systems! Therefore,
        for small systems, use simple=True, which uses a simple algorithm
        that checks all the pairs.

        """

        # casts False to empty string
        if not periodic:
            periodic = ""

        print("Computing {}...".format(struct))
        time_begin = time()

        def _bond(type1, type2, pos1, pos2):
            try:
                max_length_squared = bond_lengths_squared[type1 + ":" + type2]
            except KeyError:
                if not (type1 == "He" or type2 == "He"
                        or type1 == "Ne" or type2 == "Ne"
                        or type1 == "Ar" or type2 == "Ar"):
                    print("WARNING: pair of atom types has no bond length "
                          "defined: {} and {}".format(type1, type2))
                return False
            else:
                distance_squared = np.sum((pos2 - pos1) ** 2)
                is_bonded = distance_squared < max_length_squared
                return is_bonded

        def _angle(two_bonds: set):
            if len(two_bonds) == 4:
                return False
            elif len(two_bonds) == 3:
                return True
            elif len(two_bonds) == 2:
                print("WARNING: double bond found while computing angles!")
                return False
            else:
                raise ValueError("math broke while computing angles")

        def _dihedral(two_angles: set):
            if len(two_angles) in [5, 6]:
                return False
            elif len(two_angles) == 4:
                if any(len(a.neighbors.intersection(two_angles)) == 3 for
                       a in two_angles):
                    return False  # this is one atom bonded to all other 3
                else:
                    return True
            elif len(two_angles) == 3:
                return False  # cycle or double angle
            else:
                raise ValueError("math broke while computing dihedrals")

        xp, yp, zp = False, False, False  # periodicity
        if "x" in periodic:
            xp = True
        if "y" in periodic:
            yp = True
        if "z" in periodic:
            zp = True

        if struct == "bonds":
            self.bonds = []
            bond_lengths_squared = dict()
            for (key, value) in BOND_LENGTHS.items():
                bond_lengths_squared[key] = value ** 2
        elif struct == "angles":
            self.angles = []
        elif struct == "dihedrals":
            self.dihedrals = []

        if any(self.cell[dim][dim] > 12 for dim in [0, 1, 2]) and \
                len(self.atoms) > 1000 and not simple:

            positions = []  # this list's order is important!

            if struct == "bonds":
                for atom in self.atoms:
                    positions.append(atom.position)
            elif struct == "angles":
                for bond in self.bonds:
                    positions.append(bond.position)
            elif struct == "dihedrals":
                for angle in self.angles:
                    positions.append(angle.position)

            regions, localizers, new_pos = break_regions(positions, self.cell,
                                                         struct=struct)
            nx = len(regions)
            ny = len(regions[0])
            nz = len(regions[0][0])

            if struct == "bonds":
                for (index, atom) in enumerate(self.atoms):
                    atom.position = new_pos[index][:]
            # no need for bonds and angles to be replaced, since
            # all the atoms should be in the cell by then
            # (UNLESS the bonds weren't previously computed)

            # pairs_checked = set()  # pairs (of regions) that were checked

            for i in range(nx):
                for j in range(ny):
                    for k in range(nz):
                        for (n, index) in enumerate(regions[i][j][k]):

                            # other stuff in same region
                            if struct == "bonds":
                                this = self.atoms[index]
                                for other_index in regions[i][j][k][n + 1:]:
                                    other = self.atoms[other_index]
                                    if _bond(this.type.real, other.type.real,
                                             this.position, other.position):
                                        self.add_bond([this, other])

                            elif struct == "angles":
                                this = self.bonds[index]
                                for other_index in regions[i][j][k][n + 1:]:
                                    other = self.bonds[other_index]
                                    both_bonds = set(this.atoms + other.atoms)
                                    if _angle(both_bonds):
                                        self.add_angle(list(both_bonds))

                            elif struct == "dihedrals":
                                this = self.angles[index]
                                for other_index in regions[i][j][k][n + 1:]:
                                    other = self.angles[other_index]
                                    both_angles = set(this.atoms + other.atoms)
                                    if _dihedral(both_angles):
                                        self.add_dihedral(list(both_angles))

                            # stuff in the "next" regions
                            next_regions = [(i + 1, j - 1, k - 1),
                                            (i + 1, j - 1, k),
                                            (i + 1, j - 1, k + 1),
                                            (i + 1, j, k - 1),
                                            (i + 1, j, k),
                                            (i + 1, j, k + 1),
                                            (i + 1, j + 1, k - 1),
                                            (i + 1, j + 1, k),
                                            (i + 1, j + 1, k + 1),
                                            (i, j + 1, k - 1),
                                            (i, j + 1, k),
                                            (i, j + 1, k + 1),
                                            (i, j, k + 1)]

                            for next_region in next_regions:
                                ii, jj, kk = next_region

                                # periodic pairs
                                if periodic and (struct == "bonds"):
                                    translation = np.array([0.0, 0.0, 0.0])
                                    if (i == nx - 1) and (ii == nx) and xp:
                                        ii = 0
                                        translation += self.cell[0]
                                    elif (i == 0) and (ii == -1) and xp:
                                        ii = nx - 1
                                        translation -= self.cell[0]
                                    if (j == ny - 1) and (jj == ny) and yp:
                                        jj = 0
                                        translation += self.cell[1]
                                    elif (j == 0) and (jj == -1) and yp:
                                        jj = ny - 1
                                        translation -= self.cell[1]
                                    if (k == nz - 1) and (kk == nz) and zp:
                                        kk = 0
                                        translation += self.cell[2]
                                    elif (k == 0) and (kk == -1) and zp:
                                        kk = nz - 1
                                        translation -= self.cell[2]

                                    # avoids the same pair of regions twice
                                    # (needed only for systems short in a dim)
                                    # if (i == ii) and (j == jj) and (k == kk):
                                    #     continue
                                    # this_pair = "{}-{}-{} {}-{}-{}".format(
                                    #     i, j, k, ii, jj, kk)
                                    # same_pair = "{}-{}-{} {}-{}-{}".format(
                                    #     ii, jj, kk, i, j, k)
                                    # if this_pair in pairs_checked:
                                    #     continue
                                    # if same_pair in pairs_checked:
                                    #     continue
                                    # pairs_checked.add(this_pair)

                                    this = self.atoms[index]
                                    for other_index in regions[ii][jj][kk]:
                                        other = self.atoms[other_index]
                                        position2 = np.copy(other.position)
                                        position2 += translation
                                        if _bond(this.type.real,
                                                 other.type.real,
                                                 this.position,
                                                 position2):
                                            self.add_bond([this, other])

                                elif periodic:  # for angles and dihedrals
                                    if (i == nx - 1) and (ii == nx):
                                        ii = 0
                                    elif (i == 0) and (ii == -1):
                                        ii = nx - 1
                                    if (j == ny - 1) and (jj == ny):
                                        jj = 0
                                    elif (j == 0) and (jj == -1):
                                        jj = ny - 1
                                    if (k == nz - 1) and (kk == nz):
                                        kk = 0
                                    elif (k == 0) and (kk == -1):
                                        kk = nz - 1

                                    if struct == "angles":
                                        this = self.bonds[index]
                                        for other_index in regions[ii][jj][kk]:
                                            other = self.bonds[other_index]
                                            both_bonds = set(this.atoms + other.atoms)
                                            if _angle(both_bonds):
                                                self.add_angle(list(both_bonds))

                                    elif struct == "dihedrals":
                                        this = self.angles[index]
                                        for other_index in regions[ii][jj][kk]:
                                            other = self.angles[other_index]
                                            both_angles = set(this.atoms + other.atoms)
                                            if _dihedral(both_angles):
                                                self.add_dihedral(list(both_angles))

                                else:  # for non-periodic systems
                                    if (i == nx - 1) and (ii == nx):
                                        continue
                                    elif (i == 0) and (ii == -1):
                                        continue
                                    if (j == ny - 1) and (jj == ny):
                                        continue
                                    elif (j == 0) and (jj == -1):
                                        continue
                                    if (k == nz - 1) and (kk == nz):
                                        continue
                                    elif (k == 0) and (kk == -1):
                                        continue

                                    if struct == "bonds":
                                        this = self.atoms[index]
                                        for other_index in regions[ii][jj][kk]:
                                            other = self.atoms[other_index]
                                            if _bond(this.type.real,
                                                     other.type.real,
                                                     this.position,
                                                     other.position):
                                                self.add_bond([this, other])

                                    elif struct == "angles":
                                        this = self.bonds[index]
                                        for other_index in regions[ii][jj][kk]:
                                            other = self.bonds[other_index]
                                            bonds = set(this.atoms +
                                                        other.atoms)
                                            if _angle(bonds):
                                                self.add_angle(list(bonds))

                                    elif struct == "dihedrals":
                                        this = self.angles[index]
                                        for other_index in regions[ii][jj][kk]:
                                            other = self.angles[other_index]
                                            angles = set(this.atoms +
                                                         other.atoms)
                                            if _dihedral(angles):
                                                self.add_dihedral(list(angles))

        else:  # small systems, simple loop, slow
            if struct == "bonds":
                cut = 0.8  # this might need adjustment
                cut_x = np.abs(self.cell[0][0]) * cut
                cut_y = np.abs(self.cell[1][1]) * cut
                cut_z = np.abs(self.cell[2][2]) * cut
                for (i, atom) in enumerate(self.atoms):
                    for (j, other) in enumerate(self.atoms[i + 1:], i + 1):
                        distance = self.distance_vector(i + 1, j + 1)
                        if xp:
                            if distance[0] > cut_x:
                                distance -= self.cell[0]
                            elif distance[0] < -cut_x:
                                distance += self.cell[0]
                        if yp:
                            if distance[1] > cut_y:
                                distance -= self.cell[1]
                            elif distance[1] < -cut_y:
                                distance += self.cell[1]
                        if zp:
                            if distance[2] > cut_z:
                                distance -= self.cell[2]
                            elif distance[2] < -cut_z:
                                distance += self.cell[2]
                        if _bond(atom.type.real, other.type.real,
                                 np.array([0.0, 0.0, 0.0]), distance):
                            self.add_bond([atom, other])
            elif struct == "angles":
                self.compute_angles_simple()
            elif struct == "dihedrals":
                self.compute_dihedrals_simple()

        time_end = time()
        t = round(time_end - time_begin, 3)
        print(struct.capitalize() + " computed in {} seconds".format(t))

    def compute_angles_simple(self):
        self.angles = []
        for (index1, bond1) in enumerate(self.bonds):
            for (index2, bond2) in enumerate(self.bonds[index1 + 1:],
                                             index1 + 1):
                both_bonds = set(bond1.atoms + bond2.atoms)
                if len(both_bonds) == 4:
                    continue
                elif len(both_bonds) == 3:
                    self.add_angle(list(both_bonds))
                elif len(both_bonds) == 2:
                    print("WARNING: double bond found while computing angles!")
                else:
                    raise ValueError("math broke while computing angles")

    def compute_dihedrals_simple(self):
        self.dihedrals = []
        for (index1, angle1) in enumerate(self.angles):
            for (index2, angle2) in enumerate(self.angles[index1 + 1:],
                                              index1 + 1):
                both_angles = set(angle1.atoms + angle2.atoms)
                if len(both_angles) in [5, 6]:
                    continue
                elif len(both_angles) == 4:
                    if any(len(a.neighbors.intersection(both_angles)) == 3
                           for a in both_angles):
                        continue  # this is one atom bonded to all other 3
                    else:
                        self.add_dihedral(list(both_angles))
                elif len(both_angles) == 3:
                    continue  # cycle or double angle
                else:
                    raise ValueError("math broke while computing dihedrals")

    def compute_impropers(self):
        # improper parameters for molecular dynamics are somewhat arbitrary
        # checking and hand-editing the par file may be needed
        self.impropers = []
        for atom in self.atoms:
            if (atom.hybridization == "sp2") and (len(atom.neighbors) == 3):
                self.add_improper([atom] + list(atom.neighbors))

    def compute_cycles(self):
        # only goes up to 8-members cycles!
        for atom1 in self.atoms:
            for atom2 in atom1.neighbors:
                for atom3 in atom2.neighbors - {atom1}:
                    if atom1 in atom3.neighbors:  # cycle 3
                        self.add_cycle([atom1, atom2, atom3])
                    for atom4 in atom3.neighbors - {atom2}:
                        if atom1 in atom4.neighbors:  # cycle 4
                            self.add_cycle([atom1, atom2, atom3, atom4])
                        for atom5 in atom4.neighbors - {atom3}:
                            if atom1 in atom5.neighbors:  # cycle 5
                                self.add_cycle([atom1, atom2, atom3, atom4,
                                                atom5])
                            for atom6 in atom5.neighbors - {atom4}:
                                if atom1 in atom6.neighbors:  # 6 ...
                                    self.add_cycle([atom1, atom2, atom3,
                                                    atom4, atom5, atom6])
                                for atom7 in atom6.neighbors - {atom5}:
                                    if atom1 in atom7.neighbors:
                                        self.add_cycle([atom1, atom2, atom3,
                                                        atom4, atom5, atom6,
                                                        atom7])
                                    for atom8 in atom7.neighbors - {atom6}:
                                        if atom1 in atom8.neighbors:
                                            self.add_cycle([atom1, atom2,
                                                            atom3, atom4,
                                                            atom5, atom6,
                                                            atom7, atom8])

    def find_planar_cycles(self):
        for order in self.cycles:
            try:
                for cycle in order:
                    if all((atom.hybridization == "sp2" or
                            atom.hybridization == "other") for atom in cycle):
                        # the "other" is there because of Sulfur
                        cycle.is_planar = True
            except TypeError:  # order is None
                continue

    def compute_molecules(self):  # simple, not optimal
        self.molecules = []
        for atom in self.atoms:
            if atom.molecule is not None:
                continue
            else:
                molecule = {atom}.union(atom.neighbors)  # atoms in molecule
                is_growing = True
                while is_growing:
                    is_growing = False
                    for mol_atom in molecule:
                        if any(other_atom not in molecule
                               for other_atom in mol_atom.neighbors):
                            molecule = molecule.union(mol_atom.neighbors)
                            is_growing = True
                self.add_molecule(list(molecule))

    def pool_topological_types(self, re_compute_types=False):
        """
        Goes through previously computed topological structures, pooling
        topological types as well as indexing the structures.

        Parameters
        ----------
        re_compute_types : bool, optional
            If each topological structure's type shall be re-computed.
            For internal use, when the type of atoms are changed on the run
            e.g. due to classification. Standard is False.

        Notes
        -----
        This may take a while for large systems.

        """
        print("Pooling topological types...")

        if re_compute_types:
            for bond in self.bonds:
                bond.compute_type()
            for angle in self.angles:
                angle.compute_type()
            for dihedral in self.dihedrals:
                dihedral.compute_type()
            for improper in self.impropers:
                improper.compute_type()

        self.atom_types = []
        for (index, atom) in enumerate(self.atoms, 1):
            atom.index = index
            if atom.type not in self.atom_types:
                self.atom_types.append(atom.type)

        self.bond_types = []
        for (index, bond) in enumerate(self.bonds, 1):
            for molecule in self.molecules:
                if any(atom in molecule.atoms for atom in bond.atoms):
                    molecule.bonds.append(bond)
            bond.index = index
            if bond.type not in self.bond_types:
                self.bond_types.append(bond.type)

        self.angle_types = []
        for (index, angle) in enumerate(self.angles, 1):
            for molecule in self.molecules:
                if any(atom in molecule.atoms for atom in angle.atoms):
                    molecule.angles.append(angle)
            angle.index = index
            if angle.type not in self.angle_types:
                self.angle_types.append(angle.type)

        self.dihedral_types = []
        for (index, dihedral) in enumerate(self.dihedrals, 1):
            for molecule in self.molecules:
                if any(atom in molecule.atoms for atom in dihedral.atoms):
                    molecule.dihedrals.append(dihedral)
            dihedral.index = index
            if dihedral.type not in self.dihedral_types:
                self.dihedral_types.append(dihedral.type)

        self.improper_types = []
        for (index, improper) in enumerate(self.impropers, 1):
            for molecule in self.molecules:
                if any(atom in molecule.atoms for atom in improper.atoms):
                    molecule.impropers.append(improper)
            improper.index = index
            if improper.type not in self.improper_types:
                self.improper_types.append(improper.type)

        for (index, molecule) in enumerate(self.molecules, 1):
            molecule.index = index
            # no type computing for molecules

        print("Topological types pooled")

    def re_index_types(self):
        """Re-indexes topological types. For internal use, e.g. after
        computing dihedrals multiplicity."""

        for (index, atom_type) in enumerate(self.atom_types, 1):
            atom_type.index = index

        for (index, bond_type) in enumerate(self.bond_types, 1):
            bond_type.index = index

        for (index, angle_type) in enumerate(self.angle_types, 1):
            angle_type.index = index

        index = 1
        for dihedral_type in self.dihedral_types:
            if isinstance(dihedral_type.index, list):
                for i in range(len(dihedral_type.index)):
                    dihedral_type.index[i] = index
                    index += 1
            else:
                dihedral_type.index = index
                index += 1

        for (index, improper_type) in enumerate(self.improper_types, 1):
            improper_type.index = index

    def compute_hybridizations(self):

        # first look, based on neighbors
        for atom in self.atoms:

            if atom.type.real == "C":
                if len(atom.neighbors) == 4:
                    atom.hybridization = "sp3"
                elif len(atom.neighbors) == 3:
                    atom.hybridization = "sp2"
                elif len(atom.neighbors) == 2:
                    atom.hybridization = "sp"
                else:
                    atom.hybridization = "other"

            elif atom.type.real == "N":
                # nitrogen is more complicated than this, see below
                if len(atom.neighbors) == 3:
                    atom.hybridization = "sp3"
                elif len(atom.neighbors) == 2:
                    atom.hybridization = "sp2"
                elif len(atom.neighbors) == 1:
                    atom.hybridization = "sp"
                else:
                    atom.hybridization = "other"

            elif atom.type.real == "O":
                if len(atom.neighbors) == 2:
                    atom.hybridization = "sp3"
                elif len(atom.neighbors) == 1:
                    atom.hybridization = "sp2"
                else:
                    atom.hybridization = "other"

            else:
                atom.hybridization = "other"
                # H, S, P, etc

        # second look, based on possible resonance
        for atom in self.atoms:
            if atom.type.real == "N":
                if atom.hybridization == "sp3":  # because 3 neighbors
                    if any(a.hybridization == "sp2" for a in atom.neighbors):
                        atom.hybridization = "sp2"
                        # is this always true? NO
            elif atom.type.real == "O":
                if atom.hybridization == "sp3":  # because 2 neighbors
                    if any(a.hybridization == "sp2" for a in atom.neighbors):
                        atom.hybridization = "sp2"
                        # is this always true?

    def compute_bond_orders(self):
        pass  # in the future maybe

    def compute_atom_classification(self):
        """
        Classifies atoms based on their neighbors.

        Example
        -------
        The molecule HN=CH2 or HN=CHH has its atoms classified as:
        H : "imine H"
        N : "imine N"
        C : "imine CH2 C"
        H, H : "sp2 CH2 H"

        Other examples can be found in /files/charmm_proof, where
        classifications are comments.

        Notes
        -----
        The classification should be enough for classical force fields
        that require specific parameters, but was only verified to do
        a good job for the CHARMM General Force Field (CGenFF).

        """

        print("Classifying atoms...")

        def _classify_amine():
            """Classifies amine nitrogen."""
            # for internal use
            # note: variable n is taken from out of this function
            if ":C4" in n.topological_tags:
                n.classification = "NC4+ N"
            elif ":C3" in n.topological_tags:
                n.classification = "amine NRR N"
                for carbon in n.get_neighbors("C"):
                    if ":H3" in carbon.topological_tags:
                        carbon.classification = "amine NRR C"
                        for hydrogen in carbon.get_neighbors("H"):
                            hydrogen.classification = "amine NRR CH3 H"
            elif ":C2" in n.topological_tags:
                n.classification = "amine NHR N"
                for hydrogen in n.get_neighbors("H"):
                    hydrogen.classification = "amine NHR H"
                for carbon in n.get_neighbors("C"):
                    if ":H3" in carbon.topological_tags:
                        carbon.classification = "amine NHR C"
                        for hydrogen in carbon.get_neighbors("H"):
                            hydrogen.classification = "amine NHR CH3 H"
            elif ":C1" in n.topological_tags:
                n.classification = "amine NH2 N"
                for hydrogen in n.get_neighbors("H"):
                    hydrogen.classification = "amine NH2 H"
                for carbon in n.get_neighbors("C"):
                    if ":H3" in carbon.topological_tags:
                        carbon.classification = "amine NH2 C"
                        for hydrogen in carbon.get_neighbors("H"):
                            hydrogen.classification = "amine NH2 CH3 H"

        # gathers tags based on already computed info
        for atom in self.atoms:

            # only for carbons conjugated with carbons
            if atom.hybridization == "sp2":
                if sum((a.hybridization == "sp2") and (a.type.real == "C")
                       for a in atom.neighbors) > 1:
                    atom.topological_tags.add("conjugated")
                    for a in atom.neighbors:
                        # this is necessary at the edges
                        if a.hybridization == "sp2":
                            a.topological_tags.add("conjugated")

            # adds tags like ":H2", ":O1"
            counting = dict()
            for neighbor in atom.neighbors:
                neighbor_type = neighbor.type.real
                if neighbor_type in counting.keys():
                    counting[neighbor_type] += 1
                else:
                    counting[neighbor_type] = 1
            for (neighbor_type, count) in counting.items():
                atom.topological_tags.add(":{}{}".format(
                    neighbor_type, str(count)))

            for other in atom.get_neighbors("C"):
                if (atom.hybridization == "sp2") and (other.hybridization == "sp2"):
                    atom.topological_tags.add("=C")
                else:
                    atom.topological_tags.add("-C")

            for other in atom.get_neighbors("O"):
                if (atom.hybridization == "sp2") and (other.hybridization == "sp2"):
                    atom.topological_tags.add("=O")
                else:
                    atom.topological_tags.add("-O")

            for other in atom.get_neighbors("N"):
                if (atom.hybridization == "sp2") and (other.hybridization == "sp2"):
                    atom.topological_tags.add("=N")
                else:
                    atom.topological_tags.add("-N")

        # clears
        for atom in self.atoms:
            atom.classification = None

        # does the classification based on the tags gathered above
        # first round, classifications by C
        for atom in self.atoms:
            # obs: cycles are re-classified later
            if atom.classification is not None:
                continue
            if atom.type.real != "C":
                continue

            # CARBONATE
            elif ":O3" in atom.topological_tags:
                atom.classification = "CO3- C"
                for o in atom.get_neighbors("O"):
                    o.classification = "CO3- O"

            # GUANIDINE
            elif ":N3" in atom.topological_tags:
                n_sp2 = None
                for nitrogen in atom.get_neighbors("N"):
                    if nitrogen.hybridization == "sp2":
                        n_sp2 = nitrogen
                if n_sp2 is None:
                    atom.classification = "guanidinium C"
                    for n in atom.get_neighbors("N"):
                        n.classification = "guanidinium N"
                else:
                    atom.classification = "guanidine C"
                    for n in atom.get_neighbors("N"):
                        if len(n.neighbors) == 2:
                            n.classification = "guanidine =N"
                            h = n.get_neighbors("H")
                            if h:
                                h = h[0]
                                h.classification = "imine H"
                        else:
                            _classify_amine()

            # AMIDINE
            elif ":N2" in atom.topological_tags:
                n_sp2 = None
                for nitrogen in atom.get_neighbors("N"):
                    if nitrogen.hybridization == "sp2":
                        n_sp2 = nitrogen
                if n_sp2 is None:
                    atom.classification = "amidinium C"
                    for n in atom.get_neighbors("N"):
                        n.classification = "amidinium N"
                else:
                    atom.classification = "amidine C"
                    for n in atom.get_neighbors("N"):
                        if len(n.neighbors) == 2:
                            n.classification = "amidine =N"
                            h = n.get_neighbors("H")
                            if h:
                                h = h[0]
                                h.classification = "imine H"
                        else:
                            _classify_amine()

                if all(n.hybridization == "sp3"
                       for n in atom.get_neighbors("N")):
                    atom.classification = "amidinium C"
                else:
                    atom.classification = "amidine C"

            # AMIDE
            elif ("=O" in atom.topological_tags) and (":N1" in atom.topological_tags):
                atom.classification = "amide C"
                atom.get_neighbors("O")[0].classification = "amide O"
                n = atom.get_neighbors("N")[0]
                hydrogens = atom.get_neighbors("H")
                if hydrogens:
                    hydrogens[0].classification = "formamide H"
                if ":H2" in n.topological_tags:
                    n.classification = "amide NH2 N"
                elif ":H1" in n.topological_tags:
                    n.classification = "amide NHR N"
                else:
                    n.classification = "amide NRR' N"
                for h in n.get_neighbors("H"):
                    h.classification = "amide H"

            # IMINE and NITRO
            elif "=N" in atom.topological_tags:
                n = atom.get_neighbors("N")[0]
                if ":O2" in n.topological_tags:
                    atom.classification = "nitro C"
                    n.classification = "nitro N"
                    for o in n.get_neighbors("O"):
                        o.classification = "nitro O"
                else:
                    if ":H1" in atom.topological_tags:
                        atom.classification = "imine CRH C"
                    elif ":H2" in atom.topological_tags:
                        atom.classification = "imine CH2 C"
                    n = atom.get_neighbors("N")[0]
                    n.classification = "imine N"
                    h = n.get_neighbors("H")
                    if h:
                        h = h[0]
                        h.classification = "imine H"

            elif ":N1" in atom.topological_tags:

                # CYANIDE
                n = atom.get_neighbors("N")[0]
                if n.hybridization == "sp":
                    atom.classification = "cyanide C"
                    n.classification = "cyanide N"

                # AMINE
                else:
                    _classify_amine()

            elif ("=O" in atom.topological_tags) and (":O2" in atom.topological_tags):
                oxygens = atom.get_neighbors("O")
                assert len(oxygens) == 2
                o, double_o = None, None
                for oxygen in oxygens:
                    if len(oxygen.neighbors) == 2:
                        o = oxygen
                    else:
                        double_o = oxygen
                assert o is not None
                assert double_o is not None

                # ACID
                if ":H1" in o.topological_tags:
                    atom.classification = "acid C"
                    o.classification = "acid -O"
                    double_o.classification = "acid =O"
                    o.get_neighbors("H")[0].classification = "acid H"
                    h = atom.get_neighbors("H")  # bonded to the carbon
                    if h:
                        h = h[0]
                        h.classification = "formic acid H"

                # ESTER
                elif ":C2" in o.topological_tags:
                    atom.classification = "ester C"
                    o.classification = "ester -O"
                    double_o.classification = "ester =O"
                    for c in o.get_neighbors("C"):
                        if c is not atom:
                            c.classification = "ether C"
                            # should this be more specific?
                            for h in c.get_neighbors("H"):
                                h.classification = "ether C H"
                                # should this be more specific?
                    h = atom.get_neighbors("H")
                    if h:
                        h = h[0]
                        h.classification = "formic acid H"

            # CARBOXYLATES
            elif ":O2" in atom.topological_tags:
                # both oxygens are sp2
                if atom.hybridization == "sp":
                    atom.classification = "CO2 C"
                    for o in atom.get_neighbors("O"):
                        o.classification = "CO2 O"
                else:
                    atom.classification = "CO2- C"
                    for o in atom.get_neighbors("O"):
                        o.classification = "CO2- O"

            elif "=O" in atom.topological_tags:

                # ALDEHYDE
                if ":H1" in atom.topological_tags:
                    o = atom.get_neighbors("O")[0]
                    if len(o.neighbors) == 1:
                        atom.classification = "aldehyde C"
                        atom.get_neighbors("O")[0].classification = "aldehyde O"
                        atom.get_neighbors("H")[0].classification = "aldehyde H"
                    else:
                        atom.classification = "CX sp2 C"
                        atom.get_neighbors("H")[0].classification = "sp2 CHR H"

                # KETONE
                else:
                    atom.classification = "ketone C"
                    atom.get_neighbors("O")[0].classification = "ketone O"

            elif "-O" in atom.topological_tags:
                o = atom.get_neighbors("O")[0]

                # HYDROXYL
                if ":H1" in o.topological_tags:
                    atom.classification = "hydroxyl C"
                    o.classification = "hydroxyl O"
                    o.get_neighbors("H")[0].classification = "hydroxyl H"
                    # note: ionized alcohol oxygen was left out

                # ETHER
                elif ":C2" in o.topological_tags:
                    if atom.hybridization == "sp3":
                        atom.classification = "ether C"
                        o.classification = "ether O"
                        for c in o.get_neighbors("C"):
                            for h in c.get_neighbors("H"):
                                h.classification = "ether C H"
                    elif atom.hybridization == "sp2":
                        atom.classification = "CX sp2 C"
                        atom.get_neighbors("H")[0].classification = "sp2 CHR H"

            # FLUOR HALOALKANES
            elif ":F1" in atom.topological_tags:
                atom.classification = "CF C"
                for f in atom.get_neighbors("F"):
                    f.classification = "CF F"
                for h in atom.get_neighbors("H"):
                    h.classification = "CF H"
            elif ":F2" in atom.topological_tags:
                atom.classification = "CF2 C"
                for f in atom.get_neighbors("F"):
                    f.classification = "CF2 F"
                for h in atom.get_neighbors("H"):
                    h.classification = "CF2 H"
            elif ":F3" in atom.topological_tags:
                atom.classification = "CF3 C"
                for f in atom.get_neighbors("F"):
                    f.classification = "CF3 F"
            elif ":F4" in atom.topological_tags:
                atom.classification = "CF4 C"
                for f in atom.get_neighbors("F"):
                    f.classification = "CF4 F"

            # HYDROCARBONS
            if atom.classification is None:
                if atom.hybridization == "sp":
                    if ":H1" in atom.topological_tags:
                        atom.classification = "CH sp C"
                        atom.get_neighbors("H")[0].classification = "sp CH H"
                    else:
                        atom.classification = "sp C"
                        # should this be more specific?
                elif atom.hybridization == "sp2":
                    if "conjugated" in atom.topological_tags:
                        if ":H2" in atom.topological_tags:
                            atom.classification = "conjugated CH2 sp2 C"
                            for h in atom.get_neighbors("H"):
                                h.classification = "sp2 CH2 H"
                        elif ":H1" in atom.topological_tags:
                            atom.classification = "conjugated CHR sp2 C"
                            for h in atom.get_neighbors("H"):
                                h.classification = "sp2 CHR H"
                        else:
                            atom.classification = "conjugated CRR' sp2 C"
                    else:
                        if ":H2" in atom.topological_tags:
                            atom.classification = "CH2 sp2 C"
                            for h in atom.get_neighbors("H"):
                                h.classification = "sp2 CH2 H"
                        elif ":H1" in atom.topological_tags:
                            atom.classification = "CHR sp2 C"
                            for h in atom.get_neighbors("H"):
                                h.classification = "sp2 CHR H"
                        else:
                            atom.classification = "CRR' sp2 C"
                elif atom.hybridization == "sp3":
                    if ":H3" in atom.topological_tags:
                        atom.classification = "CH3 sp3 C"
                        for h in atom.get_neighbors("H"):
                            h.classification = "sp3 CH3 H"
                    elif ":H2" in atom.topological_tags:
                        atom.classification = "CH2 sp3 C"
                        for h in atom.get_neighbors("H"):
                            h.classification = "sp3 CH2 H"
                    elif ":H1" in atom.topological_tags:
                        atom.classification = "CH sp3 C"
                        for h in atom.get_neighbors("H"):
                            h.classification = "sp3 CH H"
                    elif ":C4" in atom.topological_tags:
                        atom.classification = "CC4 C"
                    elif ":H4" in atom.topological_tags:
                        atom.classification = "CH4 C"
                        for h in atom.get_neighbors("H"):
                            h.classification = "CH4 H"

                # other HALOALKANES
                # note that this is independent of the above classifications
                if ":Cl1" in atom.topological_tags:
                    for cl in atom.get_neighbors("Cl"):
                        cl.classification = "CCl Cl"
                elif ":Cl2" in atom.topological_tags:
                    for cl in atom.get_neighbors("Cl"):
                        cl.classification = "CCl2 Cl"
                elif ":Cl3" in atom.topological_tags:
                    for cl in atom.get_neighbors("Cl"):
                        cl.classification = "CCl3 Cl"
                elif ":Br1" in atom.topological_tags:
                    for br in atom.get_neighbors("Br"):
                        br.classification = "CBr Br"
                elif ":Br2" in atom.topological_tags:
                    for br in atom.get_neighbors("Br"):
                        br.classification = "CBr2 Br"
                elif ":Br3" in atom.topological_tags:
                    for br in atom.get_neighbors("Br"):
                        br.classification = "CBr3 Br"
                elif ":I1" in atom.topological_tags:
                    atom.get_neighbors("I")[0].classification = "CI I"

        # second round, classifications by non-C atoms
        for atom in self.atoms:

            if atom.classification is not None:
                continue

            if atom.type.real == "H":
                if ":H1" in atom.topological_tags:
                    atom.classification = "H2 H"
                    continue
                c = atom.get_neighbors("C")
                if c:
                    c = c[0]
                else:
                    continue
                if c.hybridization == "sp":
                    atom.classification = "sp CH H"
                elif c.hybridization == "sp2":
                    if "conjugated" in c.topological_tags:
                        if ":H2" in c.topological_tags:
                            atom.classification = "sp2 CH2 H"
                        elif ":H1" in c.topological_tags:
                            atom.classification = "sp2 CHR H"
                    else:
                        if ":H2" in c.topological_tags:
                            atom.classification = "sp2 CH2 H"
                        elif ":H1" in c.topological_tags:
                            atom.classification = "sp2 CHR H"
                elif c.hybridization == "sp3":
                    if ":H3" in c.topological_tags:
                        atom.classification = "sp3 CH3 H"
                    elif ":H2" in c.topological_tags:
                        atom.classification = "sp3 CH2 H"
                    elif ":H1" in c.topological_tags:
                        atom.classification = "sp3 CH H"

            elif atom.type.real == "N":
                # note: most of the protonated nitrogens were left out!
                if ":H4" in atom.topological_tags:
                    atom.classification = "NH4+ N"
                    for h in atom.get_neighbors("H"):
                        h.classification = "NH4+ H"
                elif ":H3" in atom.topological_tags:
                    atom.classification = "NH3 N"
                    for h in atom.get_neighbors("H"):
                        h.classification = "NH3 H"
                elif ":H2" in atom.topological_tags:
                    if ":N1" in atom.topological_tags:
                        atom.classification = "NH2N N"
                        for h in atom.get_neighbors("H"):
                            h.classification = "NH2N H"

            elif atom.type.real == "O":
                if ":H2" in atom.topological_tags:  # water
                    atom.classification = "H2O O"
                    for h in atom.get_neighbors("H"):
                        h.classification = "H2O H"

            elif atom.type.real == "S":
                if ":O4" in atom.topological_tags:
                    atom.classification = "SO4 S"
                    for o in atom.get_neighbors("O"):
                        o.classification = "SO4 =O"
                elif ":S1" in atom.topological_tags:
                    atom.classification = "CSSC S"
                elif ":H1" in atom.topological_tags:
                    atom.classification = "SH S"
                    atom.get_neighbors("H")[0].classification = "SH H"
                elif ":C2" in atom.topological_tags:
                    atom.classification = "CSC S"

            elif atom.type.real == "P":
                if ":O4" in atom.topological_tags:
                    if any(":P2" in o.topological_tags
                           for o in atom.get_neighbors("O")):
                        atom.classification = "pyrophosphate P"
                    else:
                        atom.classification = "PO4 P"
                    for o in atom.get_neighbors("O"):
                        o.classification = "PO4 =O"

            elif atom.type.real == "Al":
                if ":F4" in atom.topological_tags:
                    atom.classification = "AlF4 Al"
                    for f in atom.get_neighbors("F"):
                        f.classification = "AlF4 F"

        # third round, changes cycles
        # this is currently very specific for CGenFF
        for atom in self.atoms:

            if not atom.cycles:
                continue

            if atom.type.real == "C":

                if ("cycle 5" in atom.topological_tags and
                        "cycle 3" in atom.topological_tags):
                    atom.classification = "bridgehead C"

                elif ("cycle 5" in atom.topological_tags and
                      "cycle 6" in atom.topological_tags):
                    atom.classification = "bridge C"

                elif ("cycle 5" in atom.topological_tags and
                      "cycle 7" in atom.topological_tags):
                    atom.classification = "azulene bridge C"

                elif "cycle 7" in atom.topological_tags:
                    if any(c.is_planar and len(c) == 7 for c in atom.cycles):
                        atom.classification = "7-ring aromatic C"
                        atom.get_neighbors("H")[0].classification = "7-ring aromatic H"

                elif "cycle 3" in atom.topological_tags:
                    if len(atom.cycles) == 1:
                        if all(a.type.real == "C" for a in atom.cycles[0]):
                            atom.classification = "cyclopropyl C"

                # note: C for 4-members cycle was reserved, but we ignore it

                elif "cycle 5" in atom.topological_tags:
                    # CARE with 5-members cycles with N

                    if atom.hybridization == "sp2":
                        if "=N" in atom.topological_tags:
                            if (":N2" in atom.topological_tags or
                                    ":O1" in atom.topological_tags):  # etc?
                                atom.classification = "5-ring XC=N C"
                                # X is heteroatom
                            else:
                                atom.classification = "5-ring C=N C"
                        else:
                            atom.classification = "5-ring sp2 C"

                        # classifies hydrogens
                        if any(c.is_planar and (len(c) == 5)
                               for c in atom.cycles):
                            for h in atom.get_neighbors("H"):
                                if (":N1" in atom.topological_tags or
                                        ":N2" in atom.topological_tags or
                                        ":O1" in atom.topological_tags or
                                        ":S1" in atom.topological_tags):
                                    h.classification = "5-ring planar XC H"
                                    # X is heteroatom
                                else:
                                    h.classification = "5-ring planar H"

                    elif atom.hybridization == "sp3":
                        if "-N" in atom.topological_tags:
                            if ":H1" in atom.topological_tags:
                                atom.classification = "5-ring HC-N C"
                            elif ":H2" in atom.topological_tags:
                                atom.classification = "5-ring H2C-N C"
                        else:
                            if ":H1" in atom.topological_tags:
                                atom.classification = "5-ring CH C"
                            elif ":H2" in atom.topological_tags:
                                atom.classification = "5-ring CH2 C"
                            elif ":C4" in atom.topological_tags:
                                atom.classification = "5-ring CC4 C"

                elif "cycle 6" in atom.topological_tags:

                    # checks if cycle 6 is aromatic
                    planar_six_cycle = None
                    for c in atom.cycles:
                        if c.is_planar and (len(c) == 6):
                            planar_six_cycle = c
                    if not planar_six_cycle:
                        continue

                    # 6-RING AROMATICS
                    for n in atom.get_neighbors("N"):
                        if ("amine" in n.classification or
                            "imine" in n.classification) and (not n.cycles):
                            n.classification = "aniline N"  # conjugated
                            for h in n.get_neighbors("H"):
                                h.classification = "aniline H"
                    if "amide" in atom.classification:
                        atom.classification = "6-ring aromatic amide C"
                    elif "=O" in atom.topological_tags:
                        atom.get_neighbors("O")[0].classification = "6-ring aromatic C=O O"
                    elif ((":N2" in atom.topological_tags or
                           ":N3" in atom.topological_tags)
                          and "=N" in atom.topological_tags):
                        atom.classification = "6-ring aromatic NC=N C"
                        for h in atom.get_neighbors("H"):
                            h.classification = "6-ring planar XC H"
                    elif any(("cycle 6" in a.topological_tags) and
                             (a.cycles[0] is not atom.cycles[0])
                             for a in atom.get_neighbors("C")):
                        atom.classification = "biphenyl C"
                    elif any(("=O" in a.topological_tags) and
                             (a.type.real == "C") for
                             a in planar_six_cycle.atoms):
                        # CARE with the resonances in this case
                        atom.classification = "6-ring aromatic with C=O C"
                    elif ":N1" in atom.topological_tags:
                        atom.classification = "6-ring aromatic CN C"
                        for h in atom.get_neighbors("H"):
                            h.classification = "6-ring planar XC H"
                    elif ":F1" in atom.topological_tags:
                        atom.classification = "6-ring aromatic CF C"
                        atom.get_neighbors("F")[0].classification = "aromatic F"
                    elif ":Cl1" in atom.topological_tags:
                        atom.get_neighbors("Cl")[0].classification = "aromatic Cl"
                    elif ":Br1" in atom.topological_tags:
                        atom.get_neighbors("Br")[0].classification = "aromatic Br"
                    elif ":I1" in atom.topological_tags:
                        atom.get_neighbors("I")[0].classification = "aromatic I"
                    else:
                        atom.classification = "6-ring aromatic C"
                        for h in atom.get_neighbors("H"):
                            h.classification = "6-ring aromatic H"

            elif atom.type.real == "N":

                if ("cycle 5" in atom.topological_tags and
                        "cycle 6" in atom.topological_tags):
                    atom.classification = "bridge N"

                elif "cycle 5" in atom.topological_tags:

                    if "amide" in atom.classification:
                        if any("amide" in c.topological_tags and
                               "cycle 5" in c.topological_tags for
                               c in atom.get_neighbors("C")):
                            atom.classification = "5-ring amide N"
                            # both amide C and amide N in the ring
                        else:
                            atom.classification = "amide NRR' N"

                    elif atom.hybridization == "sp2":
                        # below, checks if cycle 5 is aromatic
                        planar_five_cycle = None
                        for c in atom.cycles:
                            if c.is_planar and (len(c) == 5):
                                planar_five_cycle = c
                        if not planar_five_cycle:
                            continue
                        # below, both are sp2
                        if len(atom.neighbors) == 3:
                            atom.classification = "5-ring planar 3-bond N"
                            for h in atom.get_neighbors("H"):
                                h.classification = "5-ring planar X H"
                        elif len(atom.neighbors) == 2:
                            atom.classification = "5-ring planar 2-bond N"

                    elif atom.hybridization == "sp3":
                        if ":H1" in atom.topological_tags:
                            atom.classification = "5-ring amine NH N"

                elif "cycle 6" in atom.topological_tags:
                    # below, checks if cycle 6 is aromatic
                    planar_six_cycle = None
                    for c in atom.cycles:
                        if c.is_planar and (len(c) == 6):
                            planar_six_cycle = c
                    if not planar_six_cycle:
                        continue
                    # below, both are sp2
                    if len(atom.neighbors) == 3:
                        atom.classification = "6-ring 3-bond N"
                    elif len(atom.neighbors) == 2:
                        if ((":N1" in atom.topological_tags) or
                                any(":N2" in a.topological_tags for a
                                    in atom.get_neighbors("C"))):
                            atom.classification = "6-ring 2-bond NCN N"
                        else:
                            atom.classification = "6-ring 2-bond N"

            elif atom.type.real == "O":

                if "cycle 5" in atom.topological_tags:
                    if all(a.hybridization == "sp2" for a in atom.neighbors):
                        atom.classification = "furan O"
                    else:
                        atom.classification = "5-ring ether O"

                elif "cycle 6" in atom.topological_tags:
                    if all(a.hybridization == "sp2" for a in atom.neighbors):
                        atom.classification = "pyran O"
                    else:
                        atom.classification = "6-ring ether O"

            elif atom.type.real == "S":
                if any((len(c) == 5) and c.is_planar for c in atom.cycles):
                    atom.classification = "thiophene S"

        # forth round, no specific classifications
        for atom in self.atoms:
            if atom.classification is not None:
                continue
            else:
                atom.classification = atom.type.real

        print("Atoms classified")

    def remove_molecule(self, molecule, pool_top_types=True):
        """
        Removes a molecule from the Atoms object's molecules list.
        This means all its atoms and topological structures are removed.

        Parameters
        ----------
        molecule : Molecule
            Molecule to be removed.
        pool_top_types : bool, optional
            For internal use e.g. if multiple molecules are being removed
            and the pool is held until the end. Standard is True.

        Raises
        ------
        KeyError
            If the molecule isn't in the Atoms's list of molecules.

        Notes
        -----
        Currently the cycles that are in the molecule aren't removed from
        the Atoms object, since cycles are only used for atom classification.

        """

        if molecule not in self.molecules:
            raise KeyError("molecule not in system")

        self.molecules.remove(molecule)
        for bond in molecule.bonds:
            self.bonds.remove(bond)
        for angle in molecule.angles:
            self.angles.remove(angle)
        for dihedral in molecule.dihedrals:
            self.dihedrals.remove(dihedral)
        for improper in molecule.impropers:
            self.impropers.remove(improper)
        for atom in molecule.atoms:
            self.atoms.remove(atom)

        if pool_top_types:
            self.pool_topological_types()

    def clear_cell(self):
        self._cell = None

    def disp_attrs(self):
        for key in self.__dict__:
            print("{}: {}".format(key, self.__dict__[key]))

    def count(self, atom_type):
        pass  # in the future maybe

    def add_atom(self, atom: Atom):
        """
        Adds an Atom object to this Atoms object's atom list.
        Includes its type in atom_types list, if it's not already in it.

        Parameters
        ----------
        atom : Atom
            Atom to be added.

        """
        if isinstance(atom, Atom):
            self.atoms.append(atom)
            if atom.type not in self.atom_types:
                self.atom_types.append(atom.type)
        else:
            raise TypeError("atom must be type Atom, got "
                            "{}".format(type(atom)))

    def add_atoms(self, atoms):
        """
        Given another Atoms object, adds each Atom in it to this Atoms object.
        Includes their types in atom_types list, if they're not already in it.

        Parameters
        ----------
        atoms : Atoms
            Atoms object whose atoms are to be added.

        """
        if isinstance(atoms, Atoms):
            for atom in atoms.atoms:
                self.add_atom(atom)
        else:
            raise TypeError("atoms must be type Atoms, got "
                            "{}".format(type(atoms)))

    def translate(self, vector):
        """Translates (moves) every atom by the given vector."""
        for atom in self.atoms:
            atom.translate(vector)

    def rotate(self, angle: float, axis: str):
        """
        Rotates all atoms around one axis going through the origin.

        Parameters
        ----------
        angle : float
            Angle (in degrees) of rotation.
        axis : str
            Axis to rotate around: 'x', 'y' or 'z'.

        """
        for atom in self.atoms:
            atom.rotate(angle, axis)

    def geometric_center(self):
        """Returns the geometric center of all atoms."""
        geometric_center = np.array([0.0, 0.0, 0.0])
        for atom in self.atoms:
            geometric_center += atom.position
        geometric_center /= len(self.atoms)
        return geometric_center

    def rotate_around_self(self, angle: float, axis: str):
        """
        Rotates all atoms around one axis going through their geometric center.

        Parameters
        ----------
        angle : float
            Angle (in degrees) of rotation.
        axis : str
            Axis to rotate around: 'x', 'y' or 'z'.

        """
        geometric_center = self.geometric_center()
        for atom in self.atoms:
            atom.position -= geometric_center
        self.rotate(angle, axis)
        for atom in self.atoms:
            atom.position += geometric_center

    def translate_to_zero(self):
        """Translates the system so its geometric center is at the origin."""
        self.translate(-1 * self.geometric_center())

    def translate_to_cell_center(self):
        """Translates (moves) the system so its geometric center is at the
        center of the cell."""
        if self.cell is None:
            raise NameError("cell not defined")
        else:
            self.translate_to_zero()
            cell_center = (self.cell[0] + self.cell[1] + self.cell[2]) / 2
            self.translate(cell_center)

    def extended(self, nx, ny, nz, topology=False):
        """
        Extends the system, copying the existing atoms as many times as needed.
        Returns new object, does NOT change self.

        Parameters
        ----------
        nx, ny, nz : int
            How many times the system shall be extended in x, y and z.
        topology : bool, optional
            If topology is to be extended as well. Standard is False.

        Returns
        -------
        extended : Atoms
            New Atoms object, extended from the original.

        Notes
        -----
        Meant for orthogonal boxes. Be careful when using triclinic systems.
        Does NOT compute topology internally, only extends the existing stuff.

        """

        this = self.copy(topology=topology)

        # in x
        original_0 = this.copy(topology=topology)
        for _ in range(nx - 1):
            original_0.translate([self.cell[0][0], 0, 0])
            original_0_copy = original_0.copy(topology=topology)
            this = this + original_0_copy

        # in y
        original_1 = this.copy(topology=topology)
        for _ in range(ny - 1):
            original_1.translate([0, self.cell[1][1], 0])
            original_1_copy = original_1.copy(topology=topology)
            this = this + original_1_copy

        # in z
        original_2 = this.copy(topology=topology)
        for _ in range(nz - 1):
            original_2.translate([0, 0, self.cell[2][2]])
            original_2_copy = original_2.copy(topology=topology)
            this = this + original_2_copy

        this.cell = [[nx * self.cell[0][0], 0, 0],
                     [0, ny * self.cell[1][1], 0],
                     [0, 0, nz * self.cell[2][2]]]

        extended = this
        return extended

    def filled(self, molecule_xyz_path, gap=4.0,
               centralize=False, topology=True, impropers=True,
               centralize_molecule=False, bulk_path=None):
        # honestly, PackMol does a much better job than this method
        """
        Fills the system's cell with a molecule, avoiding superposition.
        Returns new object, does not change self.

        Parameters
        ----------
        molecule_xyz_path : str
            Path to the molecule's xyz file.
        gap : float, optional
            Minimal distance between every atom in the system and every atom
            in the molecule, in angstroms. Standard is 4.0.
        centralize : bool, optional
            If the atoms in the system are to be translated to the cell center
            before filling up. Standard is False.
        topology : bool, optional
            If topology is to be computed for the molecule. If True, the whole
            molecule is removed if any atom in it superposes with the system
            (by less than the given gap). Standard is True.
        impropers : bool, optional
            If improper dihedrals are to be included in the topology.
            Standard is True.
        centralize_molecule : bool, optional
            If the molecule's atoms are to be centralized in its cell before
            filling the system with it. If True, its xyz file should have a
            Lattice. Standard is False.
        bulk_path : str, optional
            If given, is the path where an xyz file will be written with the
            molecule extended (before superposing atoms are removed).

        Returns
        -------
        filled : Atoms
            New Atoms object, with the original atoms plus filling molecules.

        Notes
        -----
        Meant to be used with orthogonal cells.

        """

        if self.cell is None:
            raise NameError("cell not defined")
        if centralize:
            self.translate_to_cell_center()

        from files.xyz import Xyz
        molecule_xyz = Xyz(molecule_xyz_path)
        if molecule_xyz.atoms.cell is None:
            raise NameError("cell in molecule xyz file not defined")
        if centralize_molecule:
            molecule_xyz.atoms.translate_to_cell_center()

        n_molecules_in_x = math.floor(self.cell[0][0] /
                                      molecule_xyz.atoms.cell[0][0])
        n_molecules_in_y = math.floor(self.cell[1][1] /
                                      molecule_xyz.atoms.cell[1][1])
        n_molecules_in_z = math.floor(self.cell[2][2] /
                                      molecule_xyz.atoms.cell[2][2])
        print("Extending molecule {}*{}*{} times...".format(
            n_molecules_in_x, n_molecules_in_y, n_molecules_in_z))
        if topology:
            molecule_xyz.atoms.compute_topology(impropers=impropers)
        bulk = molecule_xyz.atoms.extended(n_molecules_in_x, n_molecules_in_y,
                                           n_molecules_in_z, topology=topology,
                                           impropers=impropers)
        print("Molecule extended")

        if bulk_path is not None:
            bulk.write_xyz(bulk_path)

        n_total_bulk = len(bulk)
        gap_squared = gap ** 2
        molecules_to_be_removed = []
        atoms_to_be_removed = []  # just for optimizing
        print("Checking {}*{} pairs for contact...".format(len(bulk.atoms),
                                                           len(self.atoms)))
        for filling_atom in bulk.atoms:
            if filling_atom in atoms_to_be_removed:
                continue
            for standing_atom in self.atoms:
                distance_squared = np.sum((filling_atom.position -
                                           standing_atom.position) ** 2)
                if distance_squared <= gap_squared:
                    if filling_atom.molecule not in molecules_to_be_removed:
                        molecules_to_be_removed.append(filling_atom.molecule)
                        for atom in filling_atom.molecule.atoms:
                            atoms_to_be_removed.append(atom)
                    break
        print("All pairs checked")
        for molecule in molecules_to_be_removed:
            bulk.remove_molecule(molecule, pool_top_types=False)
        print("{}/{} atoms removed".format(n_total_bulk - len(bulk),
                                           n_total_bulk))

        filled = self + bulk
        filled.cell = self.cell
        if topology:
            filled.pool_topological_types()
        return filled

    def sort_atoms_by_type(self):  # alphabetical
        def return_atom_type(atom: Atom):
            return str(atom.type)
        self.atoms.sort(key=return_atom_type)

    def sort_atoms_by_index(self):
        def return_atom_index(atom: Atom):
            return atom.index
        self.atoms.sort(key=return_atom_index)

    def read_xyz(self, filename):
        """
        Reads atom types and atomic positions from given xyz file,
        as well as Lattice parameters, if there are any (one may
        use Ovito for editing an xyz file and giving it a lattice).
        Erases current atoms in the Atoms object.

        Parameters
        ----------
        filename : str
            Path to the xyz file.

        """
        # makes a dummy Xyz object
        from files.xyz import Xyz
        xyz = Xyz(filename)
        xyz.read_xyz()

        if self.atoms:
            print("WARNING: erasing current {} "
                  "atoms...".format(len(self.atoms)))
        print("Reading {} atoms from xyz file...".format(len(self.atoms)))
        self.atoms = xyz.atoms.atoms

        if xyz.atoms.cell is not None:
            self._cell = xyz.atoms.cell

    def write_xyz(self, filename, with_classification=False):
        """
        Writes xyz file with atom species and atomic positions,
        as well as cell parameters (Lattice) if there are any.

        Parameters
        ----------
        filename : str
            Path to the xyz file to be written.
        with_classification : bool, optional
            If atom classification is wanted as a comment after every line.
            Standard is False.

        """
        # makes a dummy Xyz object
        from files.xyz import Xyz
        xyz = Xyz()
        xyz.atoms.atoms = self.atoms
        if self.cell is not None:
            xyz.atoms.cell = self.cell
        xyz.write_xyz(filename, with_classification=with_classification)

    def distance_vector(self, n1, n2):
        """
        Gets the vector distance between two atoms, from atom n1 to atom n2,
        i.e. position2 minus position1. (Starts counting atoms at 1, not 0.)

        Parameters
        ----------
        n1, n2 : int
            Numbers of the first and second atoms, resp.

        Returns
        -------
        dist : numpy.array
            Distance vector from atom n1 to atom n2.

        """
        n1 -= 1
        n2 -= 1
        dist = self.atoms[n2].position - self.atoms[n1].position
        return dist

    def distance_scalar(self, n1, n2):
        """
        Gets the scalar distance between two atoms.
        (Starts counting atoms at 1, not 0.)

        Parameters
        ----------
        n1, n2 : int
            Numbers of the two atoms.

        Returns
        -------
        dist : float
            Scalar distance between the two atoms.

        """
        distance_vector = self.distance_vector(n1, n2)
        dist = np.sqrt(np.sum(distance_vector ** 2))
        return dist

    def clear_atoms(self):
        """Clears object from atoms and all topological structures."""
        self.atoms = []
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.impropers = []
        self.cycles = [None, None, None]
        self.molecules = []
        self.atom_types = []
        self.bond_types = []
        self.angle_types = []
        self.dihedral_types = []
        self.improper_types = []


# Below we define topological structures and their types, as follows:
#
#   Atoms {defined above}
#       TopologicalStructure
#           Bond
#           Angle
#           Dihedral
#           Improper
#           Cycle
#           Molecule
#
#   TopologicalType
#       BondType
#       AngleType
#       DihedralType
#       ImproperType
#
# These classes help the Atoms object to deal with atom typing
# and topological structures.


class TopologicalStructure(Atoms):
    """Class for topological structures, i.e. bonds, angles, etc."""

    def __init__(self, atoms: list):
        super().__init__(atoms)
        self.type = None  # this gets a TopologicalType
        self.index = None


class Bond(TopologicalStructure):
    """Class for bonds between two atoms."""

    def __init__(self, atoms: list):

        if len(atoms) != 2:
            raise ValueError("to Bond must be given 2 atoms, got "
                             "{}".format(len(atoms)))
        super().__init__(atoms)

        self.compute_type()

    @property
    def length(self):
        return self.distance_scalar(0, 1)

    @property
    def position(self):
        """Returns the position of the first atom."""
        # can't be the middle point due to possible periodicity
        return self.atoms[0].position

    def compute_type(self):
        bond_type_str = ":".join(str(atom.type) for atom in self.atoms)
        try:
            self.type = BondType.instances_dict[bond_type_str]
        except KeyError:
            self.type = BondType(bond_type_str)


class Angle(TopologicalStructure):
    """Class for angles among three atoms."""

    def __init__(self, atoms: list):

        if len(atoms) != 3:
            raise ValueError("to Angle must be given 3 atoms, got "
                             "{}".format(len(atoms)))
        super().__init__(atoms)
        self.check_atoms_order()
        self.compute_type()

    @property
    def angle(self):
        return NotImplemented

    @property
    def position(self):
        """Returns the position of the central atom."""
        return self.atoms[1].position

    def check_atoms_order(self):
        a0, a1, a2 = self.atoms[0], self.atoms[1], self.atoms[2]
        if all(a in a0.neighbors for a in [a1, a2]):
            self.atoms = [a1, a0, a2]
        elif all(a in a1.neighbors for a in [a0, a2]):
            self.atoms = [a0, a1, a2]
        elif all(a in a2.neighbors for a in [a0, a1]):
            self.atoms = [a0, a2, a1]
        else:
            raise ValueError("no central atom found in an angle")

    def compute_type(self):
        angle_type_str = ":".join(str(atom.type) for atom in self.atoms)
        try:
            self.type = AngleType.instances_dict[angle_type_str]
        except KeyError:
            self.type = AngleType(angle_type_str)


class Dihedral(TopologicalStructure):
    """Class for (proper) dihedrals among four atoms."""

    def __init__(self, atoms: list):

        if len(atoms) != 4:
            raise ValueError("to Dihedral must be given 4 atoms, got "
                             "{}".format(len(atoms)))
        super().__init__(atoms)
        self.is_small_cycle = False
        self.check_atoms_order()
        self.compute_type()

    @property
    def angle(self):
        return NotImplemented

    @property
    def position(self):
        """Returns the position of a central atom."""
        return self.atoms[1].position

    def check_atoms_order(self):
        central_atoms = []
        external_atoms = []

        for atom in self.atoms:
            if len(atom.neighbors.intersection(set(self.atoms))) == 2:
                central_atoms.append(atom)
            else:
                external_atoms.append(atom)

        if (len(central_atoms) != 2) or (len(external_atoms) != 2):
            if (all(any(len(c) == 3 for c in a.cycles)
                    for a in self.atoms) or
                    all(any(len(c) == 4 for c in a.cycles)
                        for a in self.atoms)):
                self.is_small_cycle = True  # 3 or 4-members cycles
            else:
                raise ValueError("central atoms not found while computing a"
                                 " dihedral: {}".format(self.atoms))

        if not self.is_small_cycle:
            if external_atoms[0] in central_atoms[0].neighbors:
                self.atoms = [external_atoms[0], central_atoms[0],
                              central_atoms[1], external_atoms[1]]
            elif external_atoms[1] in central_atoms[0].neighbors:
                self.atoms = [external_atoms[1], central_atoms[0],
                              central_atoms[1], external_atoms[0]]
            else:
                raise ValueError("right order not found for a dihedral")

    def compute_type(self):
        dihedral_type_str = ":".join(str(atom.type) for atom in self.atoms)
        try:
            self.type = DihedralType.instances_dict[dihedral_type_str]
        except KeyError:
            self.type = DihedralType(dihedral_type_str)


class Improper(TopologicalStructure):
    """Class for improper dihedrals among four atoms."""

    def __init__(self, atoms: list):

        if len(atoms) != 4:
            raise ValueError("to Improper must be given 4 atoms, got "
                             "{}".format(len(atoms)))
        super().__init__(atoms)

        self.check_atoms_order()
        self.compute_type()

    @property
    def angle(self):
        return NotImplemented

    @property
    def position(self):
        """Returns the position of the central atom."""
        return self[0].position

    def check_atoms_order(self):
        for atom in self.atoms:
            other_atoms = self.atoms[:]
            other_atoms.remove(atom)
            if all(other in atom.neighbors for other in other_atoms):
                self.atoms = [atom] + other_atoms
                break
            else:
                raise ValueError("no central atom found for an improper")

    def compute_type(self):
        improper_type_str = ":".join(str(atom.type) for atom in self.atoms)
        try:
            self.type = ImproperType.instances_dict[improper_type_str]
        except KeyError:
            self.type = ImproperType(improper_type_str)


class Cycle(TopologicalStructure):
    """Class for cycles made of three or more atoms."""

    def __init__(self, atoms: list):
        super().__init__(atoms)
        self.is_planar = False

    def __eq__(self, other):
        return set(self[:]) == set(other[:])


class Molecule(TopologicalStructure):
    """Class for molecules."""

    def __init__(self, atoms: list):
        super().__init__(atoms)

        # self.hill

    def compute_type(self):  # redefined for molecules only
        pass


class TopologicalType:
    """Class for topological types, e.g. atom types, bond types, etc."""

    # might become a metaclass or class factory some day

    def __init__(self, string):
        # string to be given as e.g. "C:C:H" or "C-C-H" for angle

        self.string = None
        if not isinstance(string, str):
            raise TypeError("TopologicalType object must be instantiated "
                            "with a string as argument, "
                            "got {}".format(type(string)))
        else:
            self.string = string

        self.order = None  # 1 for atom, 2 for bond, 3 for angle, ...
        self.list = []  # list of each atom's AtomType
        # thus, AtomTypes must be instantiated before other TopologicalTypes
        if (":" in self.string) or ("-" in self.string):
            # because AtomTypes don't need this
            self.list = self.get_list(self.string)
        else:
            self.list = [self.string]

        self.index = None
        self.params = None

        # computes self.order
        self.order = len(self.list)
        if self.order not in [1, 2, 3, 4]:
            print("WARNING: illegal topology found: order "
                  "{}".format(self.order))

    def __eq__(self, other):
        # just for overriding within each specific TopologicalType class
        return NotImplemented

    def __len__(self):
        return self.order

    def __iter__(self):
        return iter(self.list)

    def __getitem__(self, index_or_slice):
        return self.list[index_or_slice]

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return self.string

    @classmethod
    def types(cls):
        # just for overriding within each specific TopologicalType class
        return []

    @staticmethod
    def get_list(string):
        string = string.replace(" ", "").replace("-", ":")
        list_of_strings = string.split(":")
        list_of_atom_types = []
        for item in list_of_strings:
            try:
                atom_type = AtomType.instances_dict[item]
            except KeyError:
                raise TypeError("A higher-order TopologicalType has an "
                                "AtomType which hasn't been instantiated: "
                                "'{}'.".format(item, string))
            else:
                list_of_atom_types.append(atom_type)
        return list_of_atom_types


class AtomType(TopologicalType):
    """Class for atom types, e.g. 'H' or 'C2'."""

    instances_dict = dict()  # {string: TopologicalType}

    def __init__(self, string):
        super().__init__(string)

        self.mass = None
        self.charge = None
        self.real = string

        # lj parameters
        self.epsilon = None
        self.sigma = None
        self.epsilon14 = None
        self.sigma14 = None

        if self.order != 1:
            raise TypeError("AtomType must be order 1, got string "
                            "{}".format(self.string))

        self.add_instance(self)

        try:
            self.mass = ATOMIC_MASSES[self.string]
        except KeyError:  # happens when e.g. "C1" is given for "C"
            real_type = ""
            for char in self.string:
                if char.isalpha():
                    real_type += char
                else:
                    break
            try:
                self.mass = ATOMIC_MASSES[real_type]
                self.real = real_type
            except KeyError:  # happens when e.g. "CG321" is given for "C"

                # the following is an attempt to find out what real type it is
                # CARE: might go terribly wrong

                try:  # works for e.g. "AL" for "Al"
                    real_type = self.string[0].upper() + self.string[1].lower()
                    self.mass = ATOMIC_MASSES[real_type]
                    self.real = real_type
                except KeyError:  # works for e.g. "CG321" for "C"
                    try:
                        real_type = self.string[0].upper()
                        self.mass = ATOMIC_MASSES[real_type]
                        self.real = real_type
                    except KeyError:
                        if self.string == "Atom":
                            self.mass = 1.0
                        else:
                            print("WARNING: can't find the real type of "
                                  "{}".format(self.string))
                            self.mass = 1.0
                except IndexError:  # if type is very strange
                    self.mass = 1.0  # probably a dummy type

    def __str__(self, real_type=False):
        if real_type:
            return self.real
        else:
            return self.string

    def __eq__(self, other):
        if isinstance(other, AtomType):
            return self.string == other.string
        elif isinstance(other, str):
            return self.string == other
        else:
            return NotImplemented

    @classmethod
    def add_instance(cls, instance):
        if instance not in cls.instances_dict.values():
            cls.instances_dict[instance.string] = instance


class BondType(TopologicalType):
    """Class for bond types, e.g. 'C:C' or 'C7:H6'."""

    instances_dict = dict()  # {string: TopologicalType}

    # both "C:H" and "H:C" should point to the same value

    def __init__(self, string):
        super().__init__(string)

        if self.order != 2:
            raise TypeError("BondType must be order 2, got string "
                            "{}".format(self.string))

        self.bond_order = 1

        # force field parameters
        self.k = None
        self.r0 = None

        self.add_instance(self)

    def __eq__(self, other):
        if isinstance(other, BondType):
            return (self.list == other.list) or (self.list == other.list[::-1])
        elif isinstance(other, str):
            try:
                other = self.get_list(other)
            except TypeError:
                return False
            else:
                return (self.list == other) or (self.list == other[::-1])
        else:
            return NotImplemented

    @classmethod
    def add_instance(cls, instance):
        if instance not in cls.instances_dict.values():
            cls.instances_dict[":".join(str(typ) for typ
                                        in instance.list)] = instance
            cls.instances_dict[":".join(str(typ) for typ
                                        in instance.list[::-1])] = instance


class AngleType(TopologicalType):
    """Class for angle types, e.g. 'H:C:H' or 'C1:C8:C2'."""

    instances_dict = dict()  # {string: TopologicalType}

    # both "C:C:H" and "H:C:C" should point to the same value

    def __init__(self, string):
        super().__init__(string)

        if self.order != 3:
            raise TypeError("AngleType must be order 3, got string "
                            "{}".format(self.string))

        # force field parameters
        self.k = None
        self.theta0 = None
        self.k_ub = None
        self.r_ub = None

        self.add_instance(self)

    def __eq__(self, other):
        if isinstance(other, AngleType):
            return (self.list == other.list) or (self.list == other.list[::-1])
        elif isinstance(other, str):
            try:
                other = self.get_list(other)
            except TypeError:
                # if some of the atoms given is a not instantiated AtomType
                return False
            else:
                return (self.list == other) or (self.list == other[::-1])
        else:
            return NotImplemented

    @classmethod
    def add_instance(cls, instance):
        if instance not in cls.instances_dict.values():
            cls.instances_dict[":".join(str(typ) for typ
                                        in instance.list)] = instance
            cls.instances_dict[":".join(str(typ) for typ
                                        in instance.list[::-1])] = instance


class DihedralType(TopologicalType):
    """Class for (proper) dihedral types, e.g. 'H:C:C:C'."""

    instances_dict = dict()  # {string: TopologicalType}

    # both "C:C:C:H" and "H:C:C:C" should point to the same object

    def __init__(self, string):
        super().__init__(string)

        if self.order != 4:
            raise TypeError("DihedralType must be order 4, got string "
                            "{}".format(self.string))

        # force field parameters
        self.k = None
        self.n = None
        self.d = None

        self.add_instance(self)

    def __eq__(self, other):
        if isinstance(other, DihedralType):
            return (self.list == other.list) or (self.list == other.list[::-1])
        elif isinstance(other, str):
            try:
                other = self.get_list(other)
            except TypeError:
                # if some of the atoms given is a not instantiated AtomType
                return False
            else:
                return (self.list == other) or (self.list == other[::-1])
        else:
            return NotImplemented

    @classmethod
    def add_instance(cls, instance):
        if instance not in cls.instances_dict.values():
            cls.instances_dict[":".join(str(typ) for typ
                                        in instance.list)] = instance
            cls.instances_dict[":".join(str(typ) for typ
                                        in instance.list[::-1])] = instance


class ImproperType(TopologicalType):
    """Class for improper dihedral types, e.g. 'C:O:H:H' (central C)."""

    instances_dict = dict()  # {string: TopologicalType}

    # all "C:{H:N:O}" combinations should point to the same object

    def __init__(self, string):
        super().__init__(string)

        if self.order != 4:
            raise TypeError("ImproperType must be order 4, got string "
                            "{}".format(self.string))

        # force field parameters
        self.k = None
        self.x0 = 0

        self.add_instance(self)

    def __eq__(self, other):
        if isinstance(other, ImproperType):
            a0, a1, a2, a3 = self.list[0], self.list[1], \
                             self.list[2], self.list[3]
            return ([a0, a1, a2, a3] == other.list) or \
                   ([a0, a1, a3, a2] == other.list) or \
                   ([a0, a2, a1, a3] == other.list) or \
                   ([a0, a2, a3, a1] == other.list) or \
                   ([a0, a3, a1, a2] == other.list) or \
                   ([a0, a3, a2, a1] == other.list)
        elif isinstance(other, str):
            try:
                other = self.get_list(other)
            except TypeError:
                # if some of the atoms given is a not instantiated AtomType
                return False
            else:
                a0, a1, a2, a3 = self.list[0], self.list[1], \
                                 self.list[2], self.list[3]
                return ([a0, a1, a2, a3] == other) or \
                       ([a0, a1, a3, a2] == other) or \
                       ([a0, a2, a1, a3] == other) or \
                       ([a0, a2, a3, a1] == other) or \
                       ([a0, a3, a1, a2] == other) or \
                       ([a0, a3, a2, a1] == other)
        else:
            return NotImplemented

    @classmethod
    def add_instance(cls, instance):
        if instance not in cls.instances_dict.values():
            a0, a1, a2, a3 = instance.list[0], instance.list[1], \
                             instance.list[2], instance.list[3]
            cls.instances_dict[":".join([str(a0), str(a1), str(a2),
                                         str(a3)])] = instance
            cls.instances_dict[":".join([str(a0), str(a1), str(a3),
                                         str(a2)])] = instance
            cls.instances_dict[":".join([str(a0), str(a2), str(a1),
                                         str(a3)])] = instance
            cls.instances_dict[":".join([str(a0), str(a2), str(a3),
                                         str(a1)])] = instance
            cls.instances_dict[":".join([str(a0), str(a3), str(a1),
                                         str(a2)])] = instance
            cls.instances_dict[":".join([str(a0), str(a3), str(a2),
                                         str(a1)])] = instance
