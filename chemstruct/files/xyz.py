"""Class for xyz (atomic coordinates) files."""


import numpy as np
from files.main import File, find_between
from atoms import Atom, Atoms


class Xyz(File):
    """Class for xyz files."""

    def __init__(self, absolute_path=None):
        """Initiates an Xyz object.
        If absolute_path is not None, reads the xyz file."""

        super().__init__(absolute_path)
        self.tags.add("xyz")
        self.atoms = Atoms([])

        if self.absolute_path is not None:
            self.read_xyz()

    def read_xyz(self):
        """
        Reads info from xyz file, whose path was given as
        absolute_path when the object was instantiated.

        Raises
        ------
        ValueError
            If the xyz file has extra lines between atoms.
            If a Lattice is found that has the wrong number of parameters.
        TypeError
            If non-numeric values are found for an atom's position.

        """
        self._check_read()

        lines = iter(self.content)
        line = next(lines)
        n_atoms = 0
        while True:
            try:
                n_atoms = int(line.split()[0])
            except TypeError:
                line = next(lines)
                continue
            else:
                break

        lattice_indexes = self._find("Lattice")
        must_invent_cell = False
        if len(lattice_indexes) > 0:
            lattice_index = lattice_indexes[0]
            lattice = find_between(self.content[lattice_index], '"', '"')
            try:
                xx, xy, xz, yx, yy, yz, zx, zy, zz = tuple(lattice.split())
            except ValueError:
                raise ValueError("Lattice in xyz file is bad")
            self.atoms.cell = [[xx, xy, xz], [yx, yy, yz], [zx, zy, zz]]
        else:  # if no "Lattice" is found, atoms.cell is invented below
            lattice_index = 0
            must_invent_cell = True

        for line in self.content[lattice_index+1:]:
            if not line[0].isalpha():
                continue
            try:
                s, x, y, z = tuple(line.split())
            except ValueError:
                continue
                # raise ValueError("xyz file has a bad format, avoid extra lines")
            try:
                xyz = [float(x), float(y), float(z)]
            except TypeError:
                raise TypeError("bad positions in file: {} {} {}".format(x, y, z))

            atom = Atom(atom_type=s, position=xyz)
            self.atoms.add_atom(atom)

        if len(self.atoms) != n_atoms:
            print("WARNING: {} atoms declared in file, {} atoms found in file".format(
                n_atoms, len(self.atoms)))

        if must_invent_cell:
            gap = 10  # angstroms
            min_x = min(atom.position[0] for atom in self.atoms)
            max_x = max(atom.position[0] for atom in self.atoms)
            min_y = min(atom.position[1] for atom in self.atoms)
            max_y = max(atom.position[1] for atom in self.atoms)
            min_z = min(atom.position[2] for atom in self.atoms)
            max_z = max(atom.position[2] for atom in self.atoms)
            xx = max_x - min_x + gap
            yy = max_y - min_y + gap
            zz = max_z - min_z + gap
            self.atoms.cell = [[xx, 0, 0], [0, yy, 0], [0, 0, zz]]

    def write_xyz(self, filename, real_types=False,
                  with_classification=False):
        """
        Writes xyz file with info present in the Xyz object.

        Parameters
        ----------
        filename : str
            Path to output xyz file.
        real_types : bool, optional
            If real types should be used instead of types
            (e.g. C instead of C2). Standard is False.
        with_classification : bool, optional
            If atom classification is wanted as a comment after every line.
            Standard is False.

        """

        with open(filename, "w") as F:
            F.write(str(len(self.atoms)))
            F.write('\n')

            if self.atoms.cell is not None:
                F.write('Lattice="')
                F.write(' '.join([str(p) for p in self.atoms.cell[0]]) + ' ')
                F.write(' '.join([str(p) for p in self.atoms.cell[1]]) + ' ')
                F.write(' '.join([str(p) for p in self.atoms.cell[2]]) + '"' + '\n')
            else:
                F.write('\n')

            if not with_classification:
                for atom in self.atoms:
                    F.write(atom.__str__(real_type=real_types) + '\n')
            else:
                for atom in self.atoms:
                    F.write(atom.__str__(real_type=real_types) + '\t# ' +
                            atom.classification + '\n')

    def write_simple_cif(self, filename, struct_name, centralize=True):
        """
        Writes a simple cif file. Automatically diagonalizes the cell.
        Works for 1D, 2D and 3D structures, including triclinic cells.

        Parameters
        ----------
        filename : str
            Path to output cif file.
        struct_name : str
            Name of the structure, to be included in the file info.
        centralize : bool, optional
            If the atomic positions should be centralized in the cell.
            Standard is True.

        Notes
        -----
        This is a simple cif writing for theoretical (e.g. DFT) structures;
        it's not meant to be used for complete crystallographic information.

        """

        self._check_read()
        if not self.atoms.atoms:
            self.read_xyz()

        struct_name = struct_name.replace(" ", "-")

        if centralize:
            self.atoms.translate_to_cell_center()

        # the following includes the diagonalize_cell(),
        # which has to be done so the atomic coordinates correspond to the
        # xyz axes defined by the abc vectors of the cell,
        # i.e. a is parallel to x, and b is in the xy plane
        self.atoms.compute_fractional_positions()

        with open(filename, "w") as F:

            F.write("data_" + struct_name + "\n\n")

            a = np.linalg.norm(self.atoms.a)
            b = np.linalg.norm(self.atoms.b)
            c = np.linalg.norm(self.atoms.c)
            alpha = self.atoms.alpha * 180 / np.pi  # degrees
            beta = self.atoms.beta * 180 / np.pi  # degrees
            gama = self.atoms.gama * 180 / np.pi  # degrees

            F.write("_cell_length_a " + str(round(a, 4)) + "(0)" + "\n")
            F.write("_cell_length_b " + str(round(b, 4)) + "(0)" + "\n")
            F.write("_cell_length_c " + str(round(c, 4)) + "(0)" + "\n")
            F.write("_cell_angle_alpha " + str(round(alpha, 4)) + "(0)" + "\n")
            F.write("_cell_angle_beta " + str(round(beta, 4)) + "(0)" + "\n")
            F.write("_cell_angle_gamma " + str(round(gama, 4)) + "(0)" + "\n\n")

            # F.write("_symmetry_space_group_name_H-M 'P 1'" + "\n")
            # F.write("_symmetry_Int_Tables_number 1" + "\n")
            # F.write("_symmetry_cell_setting triclinic" + "\n\n")

            F.write("loop_" + "\n")
            F.write("_atom_site_label" + "\n")
            F.write("_atom_site_type_symbol" + "\n")
            F.write("_atom_site_fract_x" + "\n")
            F.write("_atom_site_fract_y" + "\n")
            F.write("_atom_site_fract_z" + "\n")

            counter = dict()
            for atom in self.atoms:

                atom_type = str(atom.type)
                if atom_type in counter.keys():
                    counter[atom_type] += 1
                else:
                    counter[atom_type] = 1

                F.write(atom_type + str(counter[atom_type]) + " " +
                        atom.type.real + " " +
                        str(round(atom.fractional_position[0], 4)) + " " +
                        str(round(atom.fractional_position[1], 4)) + " " +
                        str(round(atom.fractional_position[2], 4)) + "\n")
