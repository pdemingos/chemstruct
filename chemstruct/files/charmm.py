# Copyright 2019 Pedro G. Demingos

"""
Defines the CharmmGeneral class for dealing with
CHARMM General Force Field parameters and atomic systems.

Ref: K. Vanommeslaeghe et al., J. Comp. Chem., 2009.

"""

from files.main import File
from files.xyz import Xyz
from atoms import AtomType, BondType, AngleType, DihedralType, ImproperType
from glob import glob
from constants import PACKAGE_DIR


FROM_CLASSIFICATION_TO_CGEN = {

    # CGenFF (Charmm General Force Field)
    # X is heteroatom

    "H2O H": "HGTIP3",
    "acid H": "HGP1",
    "formic acid H": "HGR52",
    "aldehyde H": "HGR52",
    "hydroxyl H": "HGP1",
    # "ether C H": "",
    "sp CH H": "HGPAM1",  # from the site
    "sp2 CH2 H": "HGA5",
    "sp2 CHR H": "HGA4",
    "sp3 CH H": "HGA1",
    "sp3 CH2 H": "HGA2",
    "sp3 CH3 H": "HGA3",
    "NH3 H": "HGPAM3",
    "NH2N H": "HGP1",
    "amide H": "HGP1",
    "NC4+ H": "HGP5",
    "amine NH2 CH3 H": "HGAAM2",
    "amine NH2 H": "HGPAM2",
    "amine NHR CH3 H": "HGAAM1",
    "amine NHR H": "HGPAM1",
    "amine NRR CH3 H": "HGAAM0",
    "imine H": "HGP1",
    "CF H": "HGA6",
    "CF2 H": "HGA7",
    "CF3 H": "HGA7",  # approximate
    "CH4 H": "HGA3",  # approximate
    "NH4+ H": "HGP2",  # ?
    "formamide H": "HGR52",
    # "H2 H": "",
    "SH H": "HGP3",
    
    "5-ring planar XC H": "HGR52",
    "5-ring planar X H": "HGR52",
    "5-ring planar H": "HGR51",
    "6-ring planar XC H": "HGR62",
    "6-ring aromatic H": "HGR61",
    "7-ring aromatic H": "HGR71",
    "aniline H": "HGP4",

    "imine CRH C": "CG2D1",
    "imine CH2 C": "CG2D2",
    "amide C": "CG2O1",
    "acid C": "CG2O2",
    "ester C": "CG2O2",
    # "ether C"  # must be excluded for re-classification
    "aldehyde C": "CG2O4",
    "ketone C": "CG2O5",
    # "hydroxyl C"  # must be excluded for re-classification
    "CH sp C": "CG1T1",
    "cyanide C": "CG1N1",
    "sp C": "CG1T1",
    "CH2 sp2 C": "CG2D2",
    "CHR sp2 C": "CG2D1",
    "CRR' sp2 C": "CG2D1",
    "conjugated CH2 sp2 C": "CG2DC3",
    "conjugated CHR sp2 C": "CG2DC2",
    "conjugated CRR' sp2 C": "CG2DC1",
    "CH sp3 C": "CG311",
    "CH2 sp3 C": "CG321",
    "CH3 sp3 C": "CG331",
    "CH4 C": "CG331",  # approximate
    "CC4 C": "CG301",
    "CO3- C": "CG2O6",
    "guanidinium C": "CG2N1",
    "guanidine C": "CG2N1",
    "amidinium C": "CG2N2",
    "amidine C": "CG2N2",
    # "amine C": "",  # must be excluded for re-classification
    "amine NH2 C": "CG3AM2",
    "amine NHR C": "CG3AM1",
    "amine NRR C": "CG3AM0",
    "CO2 C": "CG2O7",
    "CO2- C": "CG2O3",
    "CF C": "CG322",
    "CF2 C": "CG312",
    "CF3 C": "CG302",
    "CF4 C": "CG302",  # approximate
    "CHN+ C": "",
    "CH2N+ C": "",
    "CH3N+ C": "",
    "CX sp2 C": "CG2D1O",
    
    "6-ring aromatic C": "CG2R61",
    "bridgehead C": "CG3RC1",
    "bridge C": "CG2RC0",
    "7-ring aromatic C": "CG2R71",
    "azulene bridge C": "CG2RC7",
    "cyclopropyl C": "CG3C31",
    "5-ring XC=N C": "CG2R53",
    "5-ring C=N C": "CG2R52",
    "5-ring sp2 C": "CG2R51",  # ?
    "5-ring HC-N C": "CG3C51",
    "5-ring H2C-N C": "CG3C52",
    "5-ring CH C": "CG3C51",
    "5-ring CH2 C": "CG3C52",
    "5-ring CC4 C": "CG3C50",
    "6-ring aromatic amide C": "CG2R63",
    "6-ring aromatic NC=N C": "CG2R64",
    "6-ring aromatic with C=O C": "CG2R62",  # ?
    "6-ring aromatic CN C": "CG2R61",  # ?
    "6-ring aromatic CF C": "CG2R66",
    "biphenyl C": "CG2R67",
    # no classification for cyclobutyl

    "imine N": "NG2D1",
    "amide NH2 N": "NG2S2",
    "amide NHR N": "NG2S1",
    "amide NRR' N": "NG2S0",
    "cyanide N": "NG1T1",
    "NH3 N": "NG331",
    # "NH4+ N": "",
    "NH2N N": "NG3N1",
    # "amine N": "",
    "amine NH2 N": "NG321",
    "amine NHR N": "NG311",
    "amine NRR N": "NG301",
    "guanidinium N": "NG2P1",
    "guanidine =N": "NG2D1",
    "amidinium N": "NG2P1",
    "amidine =N": "NG2D1",
    "nitro N": "NG2O1",
    "NC4+ N": "NG3P0",
    "methylamine N": "NG321",
    "dimethylamine N": "NG311",
    "trimethylamine N": "NG301",
    
    "aniline N": "NG2S3",
    "bridge N": "NG2RC0",
    "5-ring amide N": "NG2R53",
    "5-ring planar 3-bond N": "NG2R51",
    "5-ring planar 2-bond N": "NG2R50",
    "5-ring amine NH N": "NG3C51",
    "6-ring 3-bond N": "NG2R61",
    "6-ring 2-bond NCN N": "NG2R62",
    "6-ring 2-bond N": "NG2R60",

    "H2O O": "OGTIP3",
    "amide O": "OG2D1",
    "acid -O": "OG311",
    "acid =O": "OG2D1",
    "ester -O": "OG302",
    "ester =O": "OG2D1",
    "aldehyde O": "OG2D1",
    "ketone O": "OG2D3",
    "ether O": "OG301",
    "hydroxyl O": "OG311",
    "CO3- O": "OG2D2",
    "nitro O": "OG2N1",
    "CO2 O": "OG2D5",
    "CO2- O": "OG2D2",
    "SO4 -O": "OG303",
    "SO4 =O": "OG2P1",
    "PO4 -O": "OG303",
    "PO4 =O": "OG2P1",
    
    "6-ring aromatic C=O O": "OG2D4",
    "furan O": "OG2R50",
    "5-ring ether O": "OG3C51",
    "pyran O": "OG3R60",
    "6-ring ether O": "OG3C61",

    "AlF4 Al": "ALG1",
    "AlF4 F": "FGP1",
    "CF F": "FGA1",
    "CF2 F": "FGA2",
    "CF3 F": "FGA3",
    "CF4 F": "FGA3",  # approximate
    "aromatic F": "FGR1",
    "CCl Cl": "CLGA1",
    "CCl2 Cl": "CLGA1",
    "CCl3 Cl": "CLGA3",
    "aromatic Cl": "CLGR1",
    "CBr Br": "BRGA1",
    "CBr2 Br": "BRGA2",
    "CBr3 Br": "BRGA3",
    "aromatic Br": "BRGR1",
    # "CI I": "",
    "aromatic I": "IGR1",
    "SO4 S": "SG3O1",
    "CSSC S": "SG301",
    "SH S": "SG311",
    "CSC S": "SG311",
    "thiophene S": "SG2R50",
    "pyrophosphate P": "PG2",
    "PO4 P": "PG1",
    "He": "HE",
    "Ne": "NE"
}

FROM_CGENFF_TO_CHARGE = {  # suggested only!
    "HGTIP3": 0.417,
    "HGA5": 0.21,
    "HGA4": 0.15,
    "HGA1": 0.09,
    "HGA2": 0.09,
    "HGA3": 0.09,
    "CG1T1": -0.08,
    "CG1N1": 0.36,
    "CG2D2": -0.42,
    "CG2DC3": -0.42,
    "CG2DC1": -0.15,
    "CG321": -0.18,
    "CG331": -0.27,
    "CG2O7": 0.60,
    "CG2R61": -0.115,
    "NG331": -1.125,
    "OGTIP3": -0.834,
}


def charmm_proof():
    """"""
    files = glob(PACKAGE_DIR + "/files/charmm_proof/*.xyz")
    prm_file = PACKAGE_DIR + "/par_all36_cgenff.prm"
    for file in files:
        if ("proof.xyz" in file) or ("test.xyz" in file):
            continue
        print(file)
        cgenff = CharmmGeneral(xyz_file=file, par_file=prm_file)
        cgenff.compute_topology(periodic=False)
        cgenff.atoms.write_xyz(file.replace(".xyz", "_test.xyz"),
                               with_classification=True)
    for file in files:
        if ("proof.xyz" in file) or ("test.xyz" in file):
            continue
        xyz_test = Xyz(file.replace(".xyz", "_test.xyz"))
        xyz_proof = Xyz(file.replace(".xyz", "_proof.xyz"))
        for (atom_test, atom_proof) in zip(xyz_test.atoms, xyz_proof.atoms):
            if atom_test.type == atom_proof.type:
                continue
            else:
                print("BAD CLASSIFICATION FOUND BY CHARMM_PROOF: "
                      "{} in test is {} in proof, for file "
                      "'{}'".format(atom_test.type, atom_proof.type, file))
    print("CHARMM_PROOF COMPLETED")


class CharmmGeneral(File):
    """Class for a CharmmGeneral ForceField parameters file.
    Ref: K. Vanommeslaeghe et al., J. Comp. Chem., 2009.

    Reads parameters from a Charmm parameters file;
    translates the atomic classification done by an Atoms object
    to the Charmm General Force Field classification;
    writes down a par file with all relevant parameters,
    (almost) ready to be used by a LmpDat object etc."""

    def __init__(self, xyz_file=None, par_file=None):
        super().__init__()

        self.atoms = None

        if xyz_file is not None:
            self.get_xyz(xyz_file)

        self.par_file = None
        if isinstance(par_file, str):
            self.par_file = par_file

    def get_xyz(self, xyz_file):
        """
        Reads atomic positions from xyz file.

        Parameters
        ----------
        xyz_file : str
            Path to xyz file.

        """
        xyz = Xyz(xyz_file)
        self.atoms = xyz.atoms

    def get_par(self, par_file):
        """
        Stores path to Charmm General parameters file.

        Parameters
        ----------
        par_file : str
            Path to the parameters file.

        Notes
        -----
        In the CharmmGeneral FF original publication
        (K. Vanommeslaeghe et al., J. Comp. Chem., 2009)
        the parameters file is 'par_all36_cgenff.prm'.

        """
        self.par_file = par_file

    def write_xyz(self, output_path):
        xyz = Xyz()
        xyz.atoms = self.atoms
        xyz.write_xyz(output_path)

    def compute_topology(self, periodic="", simple=False):
        """
        Computes the molecular topology (bonds, angles, etc)
        for the atoms previously read (e.g. by calling get_xyz()),
        classifying atoms according to topology and the
        CHARMM General Force Field definitions.

        Parameters
        ----------
        periodic : str, optional
            Axes in which the system is periodic (e.g. 'xyz', 'xy', 'z';
            use an empty string '' for non-periodic). Standard is False.
        simple : bool, optional
            If a simple, less optimised algorithm is wanted for computing
            bonds, angles and dihedrals. (For small systems. See Notes.)
            Standard is False.

        Notes
        -----
        Computing the topology of large systems may take a while.

        """

        if self.atoms is None:
            raise NameError("No atoms")
        self.atoms.compute_topology(periodic=periodic, complete=True,
                                    hold_pool_top_types=True,
                                    simple=simple)

        # classification to types

        bad_carbons = ["hydroxyl C", "ether C", "amine C"]

        for atom in self.atoms.atoms:

            try:
                atom.type = FROM_CLASSIFICATION_TO_CGEN[atom.classification]

            # reassigns some classifications for CGenFF
            except KeyError:
                if atom.classification in bad_carbons:
                    if atom.hybridization == "sp3":
                        if ":H" in atom.topological_tags:
                            atom.classification = "CH sp3 C"
                            for h in atom.get_neighbors("H"):
                                h.classification = "sp3 CH H"
                        elif ":H2" in atom.topological_tags:
                            atom.classification = "CH2 sp3 C"
                            for h in atom.get_neighbors("H"):
                                h.classification = "sp3 CH2 H"
                        elif ":H3" in atom.topological_tags:
                            atom.classification = "CH3 sp3 C"
                            for h in atom.get_neighbors("H"):
                                h.classification = "sp3 CH3 H"
                    elif atom.hybridization == "sp2":
                        if ":H" in atom.topological_tags:
                            atom.classification = "CHR sp2 C"
                            for h in atom.get_neighbors("H"):
                                h.classification = "sp2 CHR H"

                try:
                    atom.type = FROM_CLASSIFICATION_TO_CGEN[atom.classification]
                except KeyError:
                    print("WARNING: bad classification for CGenFF: '{}' with [{}]".format(
                        atom.classification, " ".join(list(atom.topological_tags))))

        # pools
        self.atoms.pool_topological_types(re_compute_types=True)
        self.get_params()
        self.suggest_charges()

    def suggest_charges(self):
        for atom in self.atoms:
            try:
                atom.type.charge = FROM_CGENFF_TO_CHARGE[str(atom.type)]
            except KeyError:
                continue
        print("WARNING: CHECK THE SUGGESTED CHARGES IN THE PAR FILE !!!!!!!")

    def get_params(self):
        """Goes through the parameters file, assigning parameter
        values to previously computed topological structures."""
        self.absolute_path = self.par_file
        self.read_file()

        # finds keywords
        atoms_index, bonds_index, angles_index = None, None, None
        dihedrals_index, impropers_index, nonbonded_index = None, None, None
        for (index, line) in enumerate(self.content):
            if line.startswith("ATOMS"):
                atoms_index = index
            elif line.startswith("BONDS"):
                bonds_index = index
            elif line.startswith("ANGLES"):
                angles_index = index
            elif line.startswith("DIHEDRALS"):
                dihedrals_index = index
            elif line.startswith("IMPROPERS"):
                impropers_index = index
            elif line.startswith("NONBONDED"):
                nonbonded_index = index
            else:
                continue
        assert bonds_index > atoms_index
        assert angles_index > bonds_index
        assert dihedrals_index > angles_index
        assert impropers_index > dihedrals_index
        assert nonbonded_index > impropers_index

        # below, looks for Topological Types that exist in self.atoms
        # and edits their properties based on the parameters file

        # reads ATOMS parameters
        for line in self.content[atoms_index:bonds_index]:

            try:
                _, type_index, atom, mass, *comment = tuple(line.split())
            except ValueError:
                continue

            try:
                atom_type = AtomType.instances_dict[atom]
            except KeyError:
                continue

            atom_type.mass = float(mass)

        # reads BONDS parameters
        for line in self.content[bonds_index:angles_index]:

            try:
                atom1, atom2, k, r0, *comment = tuple(line.split())
            except ValueError:
                continue

            try:
                bond_type = BondType.instances_dict[atom1 + ":" + atom2]
            except KeyError:
                try:  # redundant
                    bond_type = BondType.instances_dict[atom2 + ":" + atom1]
                except KeyError:
                    continue

            bond_type.k = float(k)
            bond_type.r0 = float(r0)

        # reads ANGLES parameters
        for line in self.content[angles_index:dihedrals_index]:

            try:
                atom1, atom2, atom3, k, theta0, *rest = tuple(line.split())
            except ValueError:
                continue
            k_ub, r_ub = 0.0, 0.0
            if rest[0] != "!":
                k_ub, r_ub, *comments = rest

            try:
                angle_type = AngleType.instances_dict[atom1 + ":" + atom2 + ":" + atom3]
            except KeyError:
                try:  # redundant
                    angle_type = AngleType.instances_dict[atom3 + ":" + atom2 + ":" + atom1]
                except KeyError:
                    continue

            angle_type.k = float(k)
            angle_type.theta0 = float(theta0)
            angle_type.k_ub = float(k_ub)
            angle_type.r_ub = float(r_ub)

        # reads DIHEDRALS parameters
        for line in self.content[dihedrals_index:impropers_index]:

            try:
                atom1, atom2, atom3, atom4, k, n, d, *comments = tuple(line.split())
            except ValueError:
                continue

            try:
                dihedral_type = DihedralType.instances_dict[atom1 + ":" + atom2 + ":" + atom3 + ":" + atom4]
            except KeyError:
                try:  # redundant
                    dihedral_type = DihedralType.instances_dict[atom4 + ":" + atom3 + ":" + atom2 + ":" + atom1]
                except KeyError:
                    continue

            try:
                dihedral_type.k.append(float(k))
                dihedral_type.n.append(int(n))
                dihedral_type.d.append(int(float(d)))
            except AttributeError:
                dihedral_type.k = [float(k)]
                dihedral_type.n = [int(n)]
                dihedral_type.d = [int(float(d))]

        # reads IMPROPERS parameters
        for line in self.content[impropers_index:nonbonded_index]:

            try:
                atom1, atom2, atom3, atom4, k, *rest = tuple(line.split())
            except ValueError:
                continue

            try:
                improper_type = ImproperType.instances_dict[atom1 + ":" + atom2 + ":" + atom3 + ":" + atom4]
            except KeyError:
                continue  # nothing redundant included for impropers

            improper_type.k = float(k)
            # improper_type.x0 = 0  # it should be zero

        # reads NONBONDED parameters
        r_min_to_sigma = (1 / 2) ** (1 / 6)
        for line in self.content[nonbonded_index + 1:]:

            if line.startswith("!"):
                continue  # this is a commentary
            if line.startswith("HBOND CUTHB"):
                break  # this is the end
            if line.startswith("cutnb"):
                continue

            try:
                atom, check_zero, minus_epsilon, r_min_over2, *rest = tuple(line.split())
            except ValueError:
                continue
            if float(check_zero) != 0.0:
                # this key be 0.0 means the next term is -epsilon
                print("WARNING: 0.0 key not found in CHARMM NONBONDED parameters!")
            minus_epsilon14, r14_min_over2 = minus_epsilon, r_min_over2
            if rest[0] != "!":
                try:
                    _, minus_epsilon14, r14_min_over2, *comment = rest
                except ValueError:
                    print(rest)

            try:
                atom_type = AtomType.instances_dict[atom]
            except KeyError:
                continue

            atom_type.epsilon = round(float(minus_epsilon.replace("-", "")), 5)
            atom_type.sigma = round(float(r_min_over2) * 2 * r_min_to_sigma, 5)
            atom_type.epsilon14 = round(float(minus_epsilon14.replace("-", "")), 5)
            atom_type.sigma14 = round(float(r14_min_over2) * 2 * r_min_to_sigma, 5)

    def write_par(self, output_path):
        """
        Writes a parameters (.par) file with the topological structures
        found in the atomic system, as well as their parameters.

        Parameters
        ----------
        output_path : str
            Path to output parameters file.

        Notes
        -----
        If no parameters are found for a topological structure, they are
        left as 'None'. Then, manually adding parameters may be necessary.

        Example
        -------
        >> from files.charmm import CharmmGeneral
        >> xyz_file = "/home/mydir/pvp3.xyz"
        >> par_file = "/home/mydir/par_all36_cgenff.prm"
        >> charmm = CharmmGeneral(xyz_file=xyz_file, par_file=par_file)
        >> charmm.compute_topology()  # may take a while
        >> charmm.write_par("/home/mydir/pvp3.par")

        """

        with open(output_path, "w") as F:

            F.write("Atom Types\n\n")
            F.write("# CARE: THE CHARGES BELOW ARE SUGGESTIONS ONLY,\n")
            F.write("# CHECK IF THEY ARE REASONABLE (e.g. NEUTRAL SYSTEM)\n")
            F.write("# type\tcharge\tepsilon\tsigma\tepsil14\tsigma14\n")
            for atom_type in self.atoms.atom_types:
                F.write("\t".join([str(atom_type), str(atom_type.charge),
                                   str(atom_type.epsilon), str(atom_type.sigma),
                                   str(atom_type.epsilon14), str(atom_type.sigma14)]) + "\n")

            F.write("\nBond Types\n\n")
            F.write("# bond_type\tk\tr0\n")
            for bond_type in self.atoms.bond_types:
                F.write("\t".join([str(bond_type), str(bond_type.k), str(bond_type.r0)]) + "\n")

            F.write("\nAngle Types\n\n")
            F.write("# angle_type\t\tk\ttheta0\tk_ub\tr_ub\n")
            for angle_type in self.atoms.angle_types:
                F.write("\t".join([str(angle_type),
                                   str(angle_type.k), str(angle_type.theta0),
                                   str(angle_type.k_ub), str(angle_type.r_ub)]) + "\n")

            F.write("\nDihedral Types\n\n")
            F.write("# dihedral_type\t\tk\tn\td\tweight=1.0\n")
            for dihedral_type in self.atoms.dihedral_types:
                try:
                    for i in range(len(dihedral_type.k)):
                        F.write("\t".join([str(dihedral_type), str(dihedral_type.k[i]),
                                           str(dihedral_type.n[i]), str(dihedral_type.d[i]),
                                           "1.0"]) + "\n")
                except TypeError:  # if this dihedral type has no parameters
                    F.write(str(dihedral_type) + "\tNone\tNone\tNone\t1.0" + "\n")

            F.write("\nImproper Types\n\n")
            F.write("# improper_type\t\tk\tx0=0\n")
            for improper_type in self.atoms.improper_types:
                F.write("\t".join([str(improper_type), str(improper_type.k),
                                   str(improper_type.x0)]) + "\n")
