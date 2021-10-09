"""
Defines the CharmmGeneral class for dealing with
CHARMM General Force Field parameters and atomic systems.

Ref: K. Vanommeslaeghe et al., J. Comp. Chem., 2009.

"""

import re

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
    # "NR4+ H": "",  # ? should be more specific?
    # "NH H": "",

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
    "conjugated sp C": "CG1T1",
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
    # "CHN+ C": "",
    # "CH2N+ C": "",
    # "CH3N+ C": "",
    "CX sp2 C": "CG2D1O",
    "nitro C": "",  # ?

    "6-ring aromatic C": "CG2R61",
    "bridgehead C": "CG3RC1",
    "bridge C": "CG2RC0",
    "7-ring aromatic C": "CG2R71",
    "azulene bridge C": "CG2RC7",
    "3-ring sp3 C": "CG3C31",
    "4-ring sp3 C": "CG3C31",
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
    "NR4+ N": "NG3P0",
    "methylamine N": "NG321",
    "dimethylamine N": "NG311",
    "trimethylamine N": "NG301",
    # "two-neighbor N": "",

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
    "aromatic S": "SG2R50",  # thiophene S
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

FROM_CLASSIFICATION_TO_GAFF = {

    # GAFF (General Amber Force Field)
    # X is heteroatom
    # electrwd. groups are: N, O, F, Cl, Br

    "H2O H": "hw",
    "acid H": "ho",
    "formic acid H": "h2",
    "aldehyde H": "ha",  # ? example from the paper
    "hydroxyl H": "ho",
    "ether C H": "h1",
    # "sp CH H": "",  # reclassification
    # "sp2 CH2 H": "",  # reclassification
    # "sp2 CHR H": "",  # reclassification
    # "sp3 CH H": "",  # reclassification
    # "sp3 CH2 H": "",  # reclassification
    # "sp3 CH3 H": "",  # reclassification
    "NH3 H": "hn",
    "NH2N H": "hn",
    "amide H": "hn",
    "NC4+ H": "hn",
    "amine NH2 CH3 H": "h1",
    "amine NH2 H": "hn",
    "amine NHR CH3 H": "h1",
    "amine NHR H": "hn",
    "amine NRR CH3 H": "h1",
    "imine H": "hn",
    "CF H": "h1",
    "CF2 H": "h2",
    "CF3 H": "h3",
    "CH4 H": "hc",
    "NH4+ H": "hn",
    "formamide H": "h2",
    # "H2 H": "",
    "SH H": "hs",
    "PH H": "hp",
    "NR4+ H": "hn",
    "NH H": "hn",

    "5-ring planar XC H": "ha",  # ?
    "5-ring planar X H": "hn",  # h bonded to nitrogen
    "5-ring planar H": "ha",
    "6-ring planar XC H": "ha",  # ?
    "6-ring aromatic H": "ha",
    "7-ring aromatic H": "ha",
    "aniline H": "hn",

    "imine CRH C": "c2",
    "imine CH2 C": "c2",
    "amide C": "c",
    "acid C": "c",
    "ester C": "c",
    # "ether C": "c3",  # reclassification
    "aldehyde C": "c",
    "ketone C": "c",
    # "hydroxyl C": "c3",  # reclassification
    "CH sp C": "c1",
    "cyanide C": "c1",
    # "sp C": "",  # reclassification
    # "conjugated sp C": "",  # reclassification
    "CH2 sp2 C": "c2",
    "CHR sp2 C": "c2",
    "CRR' sp2 C": "c2",
    "conjugated CH2 sp2 C": "c2",
    # "conjugated CHR sp2 C": "",  # reclassification
    # "conjugated CRR' sp2 C": "",  # reclassification
    "CH sp3 C": "c3",
    "CH2 sp3 C": "c3",
    "CH3 sp3 C": "c3",
    "CH4 C": "c3",
    "CC4 C": "c3",
    # "CO3- C": "",
    # "guanidinium C": "",  # ?
    "guanidine C": "cz",
    # "amidinium C": "",
    # "amidine C": "",
    # "amine C": "c3",  # reclassification
    # "amine NH2 C": "c3",  # reclassification
    # "amine NHR C": "c3",  # reclassification
    # "amine NRR C": "c3",  # reclassification
    # "CO2 C": "",
    # "CO2- C": "",
    # "CF C": "",  # reclassification
    # "CF2 C": "",  # reclassification
    # "CF3 C": "",  # reclassification
    # "CF4 C": "",  # reclassification
    # "CHN+ C": "",
    # "CH2N+ C": "",
    # "CH3N+ C": "",
    # "CX sp2 C": "cs",  # ???
    "C=S C": "c",  # example from the paper
    # "nitro C": "",  # reclassification

    # "6-ring aromatic C": "",  # reclassification
    # "bridgehead C": "",  # reclassification
    # "bridge C": "",  # reclassification
    "7-ring aromatic C": "ca",  # approximate
    # "azulene bridge C": "",  # reclassification
    "3-ring sp3 C": "cx",
    "3-ring sp2 C": "cu",
    "4-ring sp3 C": "cy",
    "4-ring sp2 C": "cv",
    # "5-ring XC=N C": "",  # reclassification
    # "5-ring C=N C": "",  # reclassification
    # "5-ring sp2 C": "",  # reclassification
    # "5-ring HC-N C": "",  # reclassification
    # "5-ring H2C-N C": "",  # reclassification
    # "5-ring CH C": "",  # reclassification
    # "5-ring CH2 C": "",  # reclassification
    # "5-ring CC4 C": "",  # reclassification
    "6-ring aromatic amide C": "c",
    # "6-ring aromatic NC=N C": "",  # reclassification
    # "6-ring aromatic with C=O C": "",  # reclassification
    # "6-ring aromatic CN C": "",  # reclassification
    # "6-ring aromatic CF C": "",  # reclassification
    # "biphenyl C": "",  # reclassification

    # "imine N": "",
    "amide NH2 N": "n",
    "amide NHR N": "n",
    "amide NRR' N": "n",
    "cyanide N": "n1",
    "NH3 N": "n9",
    "NH4+ N": "n+",
    "NH2N N": "n3",
    # "amine N": "",  # reclassification
    # "amine NH2 N": "",  # reclassification
    # "amine NHR N": "",  # reclassification
    # "amine NRR N": "",  # reclassification
    # "guanidinium N": "",
    "guanidine =N": "n2",
    # "amidinium N": "",
    "amidine =N": "n2",
    "nitro N": "no",
    "NR4+ N": "n4",
    # "methylamine N": "n3",  # reclassification
    # "dimethylamine N": "n3",  # reclassification
    # "trimethylamine N": "n3",  # reclassification
    "two-neighbor N": "n2",

    "aniline N": "nh",
    # "bridge N": "",  # ?
    "5-ring amide N": "n",
    "5-ring planar 3-bond N": "na",
    # "5-ring planar 2-bond N": "",  # reclassification
    "5-ring amine NH N": "nh",
    # "6-ring 3-bond N": "",  # reclassification
    # "6-ring 2-bond NCN N": "",  # ?
    # "6-ring 2-bond N": "",  # ?

    "H2O O": "ow",
    "amide O": "o",
    "acid -O": "oh",  # hydroxyl
    "acid =O": "o",
    "ester -O": "os",
    "ester =O": "o",
    "aldehyde O": "o",
    "ketone O": "o",
    "ether O": "os",
    "hydroxyl O": "oh",
    "CO3- O": "o",  # ?
    "nitro O": "o",
    "CO2 O": "o",  # ?
    "CO2- O": "o",  # ?
    "SO4 -O": "o",  # ?
    "SO4 =O": "o",
    "PO4 -O": "o",  # ?
    "PO4 =O": "o",
    "SO O": "o",  # example from the paper
    "PO O": "o",  # example from the paper
    "POC O": "os",  # example from the paper

    "6-ring aromatic C=O O": "o",
    "furan O": "os",
    "5-ring ether O": "os",
    "pyran O": "os",
    "6-ring ether O": "os",

    "pyrophosphate P": "py",  # ? p5 in conjugated systems
    "PO4 P": "p5",  # four connected atoms
    "two-neighbor P": "p2",
    "three-neighbor conjugated P=O P": "px",
    "three-neighbor P=O P": "p4",
    "three-neighbor P": "p3",
    "four-neighbor conjugated P": "py",
    "four-neighbor P": "p5",
    # "conjugated P": "pe",  # reclassification
    # "aromatic P": "pd",  # reclassification

    "SO4 S": "s6",  # four connected atoms
    "CSSC S": "ss",  # example from the paper
    "SH S": "sh",  # ? sp3
    "CSC S": "ss",  # thio-ether
    "C=S S": "s2",  # example from the paper
    "aromatic S": "ss",  # thio-ether
    "one-neighbor S": "s",
    "two-neighbor S": "s2",
    "three-neighbor conjugated S": "sx",
    "three-neighbor S": "s4",
    "four-neighbor conjugated S": "sy",
    "four-neighbor S": "s6",

    # "AlF4 Al": "",
    "AlF4 F": "f",
    "CF F": "f",
    "CF2 F": "f",
    "CF3 F": "f",
    "CF4 F": "f",
    "aromatic F": "f",
    "CCl Cl": "cl",
    "CCl2 Cl": "cl",
    "CCl3 Cl": "cl",
    "aromatic Cl": "cl",
    "CBr Br": "br",
    "CBr2 Br": "br",
    "CBr3 Br": "br",
    "aromatic Br": "br",
    "CI I": "i",
    "aromatic I": "i",
    # "He": "",
    # "Ne": ""
}

FROM_CLASSIFICATION_TO_COMPASS = {

    # PDMS
    "sp3 CH3 H": "h1",
    "CH3 sp3 C": "c4",
    "SiO O": "o2",
    "Si": "si4c",

    # GO
    "conjugated CHR sp2 C": "c3a",
    "conjugated CRR' sp2 C": "c3a",
    "6-ring aromatic C": "c3a",
    "hydroxyl C": "c44",
    "ether C": "c44",
    "6-ring aromatic H": "h1",
    "sp2 CHR H": "h1",
    "hydroxyl H": "ho2",
    "ether O": "o_2",
    "hydroxyl O": "o_2"

}


def wildcard(type_dict, wildcard_string, wild_char="X"):
    """"""

    types = []
    pattern = re.compile(wildcard_string.replace(wild_char, "."))

    for (string, typ) in type_dict.items():
        if bool(pattern.search(string)):

            if not typ.k:  # ignores dihedrals that already have parameters
                types.append(typ)

    return types


def charmm_proof():
    """"""  # TODO write docstring
    files = glob(PACKAGE_DIR + "/files/charmm_proof/*.xyz")
    for file in files:
        if ("proof.xyz" in file) or ("test.xyz" in file):
            continue
        print(file)
        cgenff = CharmmGeneral(xyz_file=file)
        cgenff.compute_topology(periodic="")
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


def amber_proof():
    """"""  # TODO write docstring
    files = glob(PACKAGE_DIR + "/files/amber_proof/*.xyz")
    for file in files:
        if ("proof.xyz" in file) or ("test.xyz" in file):
            continue
        print(file)
        gaff = GeneralAmber(xyz_file=file)
        gaff.compute_topology(periodic="")
        gaff.atoms.write_xyz(file.replace(".xyz", "_test.xyz"),
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
    print("AMBER_PROOF COMPLETED")


class ForceField(File):
    """"""

    def __init__(self, xyz_file=None):
        super().__init__()

        self.atoms = None
        self.par_file = None
        if xyz_file is not None:
            self.get_xyz(xyz_file)

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

    def write_xyz(self, output_path):
        xyz = Xyz()
        xyz.atoms = self.atoms
        xyz.write_xyz(output_path)

    def check_water_angles(self, model="tip3p"):
        """TIP3P etc"""  # TODO fix and write docstring

        if model == "tip3p":
            water_angle = 104.52  # degrees
            water_bond = 0.9572  # angstrom
        else:
            raise ValueError("Water model {} not implemented".format(model))

        for angle in self.atoms.angles:
            if str(angle.type) == "HGTIP3:OGTIP3:HGTIP3":
                angle.set_angle(water_angle)
        for bond in self.atoms.bonds:
            if str(bond.type) in ["HGTIP3:OGTIP3", "OGTIP3:HGTIP3"]:
                bond.set_length(water_bond)

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
        >> from files.ff import CharmmGeneral
        >> xyz_file = "/home/mydir/pvp3.xyz"
        >> par_file = "/home/mydir/par_all36_cgenff.prm"
        >> charmm = CharmmGeneral(xyz_file=xyz_file, par_file=par_file)
        >> charmm.compute_topology()  # may take a while
        >> charmm.write_par("/home/mydir/pvp3.par")

        """

        with open(output_path, "w") as F:

            F.write("Atom Types\n\n")
            if "cgenff" in self.tags:
                F.write("# CARE: THE CHARGES BELOW ARE SUGGESTIONS ONLY,\n")
                F.write("# CHECK IF THEY ARE REASONABLE (e.g. NEUTRAL SYSTEM)\n")
                F.write("# type\tcharge\tepsilon\tsigma\tepsil14\tsigma14\n")
                for atom_type in self.atoms.atom_types:
                    if atom_type.epsilon14 is not None:
                        F.write("\t".join([str(atom_type), str(atom_type.charge),
                                           str(atom_type.epsilon), str(atom_type.sigma),
                                           str(atom_type.epsilon14), str(atom_type.sigma14)]) + "\n")
            else:
                F.write("# type\tcharge\tepsilon\tsigma\n")
                for atom_type in self.atoms.atom_types:
                    F.write("\t".join([str(atom_type), str(atom_type.charge),
                                       str(atom_type.epsilon), str(atom_type.sigma)]) + "\n")

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
            if "cgenff" in self.tags:  # charmm
                F.write("# improper_type\t\tk\tx0=0\n")
                for improper_type in self.atoms.improper_types:
                    F.write("\t".join([str(improper_type), str(improper_type.k),
                                       str(improper_type.x0)]) + "\n")
            else:  # amber
                F.write("# improper_type\t\tk\td\tn\n")
                for improper_type in self.atoms.improper_types:
                    F.write("\t".join([str(improper_type), str(improper_type.k),
                                       str(improper_type.d), str(improper_type.n)]) + "\n")


class CharmmGeneral(ForceField):
    """Class for a Charmm General Force Field parameters file.
    Ref: K. Vanommeslaeghe et al., J. Comp. Chem., 2009.

    Reads parameters from a Charmm parameters file;
    translates the atomic classification done by an Atoms object
    to the Charmm General Force Field classification;
    writes down a par file with all relevant parameters,
    (almost) ready to be used by a LmpDat object etc.

    In the CharmmGeneral FF original publication
    (K. Vanommeslaeghe et al., J. Comp. Chem., 2009)
    the parameters file is 'par_all36_cgenff.prm'."""

    def __init__(self, xyz_file=None):
        super().__init__(xyz_file=xyz_file)
        self.tags.add("cgenff")
        self.par_file = PACKAGE_DIR + "/par_all36_cgenff.prm"

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

        self.check_water_angles()  # bugs if molecules are periodic!
        # TODO fix this water thing

    def suggest_charges(self):
        for atom in self.atoms:
            try:
                atom.type.charge = FROM_CGENFF_TO_CHARGE[str(atom.type)]
            except KeyError:
                continue
        print("WARNING: CHECK THE SUGGESTED CHARGES IN THE PAR FILE !!!!!!")

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

            bond = atom1 + ":" + atom2
            try:
                bond_type = BondType.instances_dict[bond]
            except KeyError:
                try:  # redundant
                    bond_type = BondType.instances_dict[bond]
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

            angle = atom1 + ":" + atom2 + ":" + atom3
            try:
                angle_type = AngleType.instances_dict[angle]
            except KeyError:
                try:  # redundant
                    angle_type = AngleType.instances_dict[angle]
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

            dihedral = atom1 + ":" + atom2 + ":" + atom3 + ":" + atom4
            try:
                dihedral_type = DihedralType.instances_dict[dihedral]
            except KeyError:
                try:  # redundant
                    dihedral_type = DihedralType.instances_dict[dihedral]
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

            improper = atom1 + ":" + atom2 + ":" + atom3 + ":" + atom4
            try:
                improper_type = ImproperType.instances_dict[improper]
            except KeyError:
                improper_types = wildcard(ImproperType.instances_dict, improper)
                for improper_type in improper_types:
                    improper_type.k = float(k)
                    # improper_type.x0 = 0  # it should be zero
            else:
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


class GeneralAmber(ForceField):
    """"""

    def __init__(self, xyz_file=None):
        super().__init__(xyz_file=xyz_file)
        self.tags.add("gaff")
        self.par_file = PACKAGE_DIR + "/gaff2.dat"

    def compute_topology(self, periodic="", simple=False):
        """"""
        if self.atoms is None:
            raise NameError("No atoms")
        self.atoms.compute_topology(periodic=periodic, complete=True,
                                    hold_pool_top_types=True,
                                    simple=simple)

        # first round, classification to types
        for atom in self.atoms.atoms:
            try:
                atom.type = FROM_CLASSIFICATION_TO_GAFF[atom.classification]
                atom.tags.add("ff class done")
            except KeyError:
                continue

        bad_carbons = ["ether C", "hydroxyl C", "amine C", "amine NH2 C",
                       "amine NHR C", "amine NRR C", "CF C", "CF2 C", "CF3 C",
                       "CF4 C", "conjugated CHR sp2 C", "conjugated CRR' sp2 C",
                       "sp C", "nitro C", "conjugated sp C"]

        bad_aromatic_carbons = ["6-ring aromatic C", "bridgehead C", "bridge C",
                                "azulene bridge C", "5-ring XC=N C", "5-ring C=N C",
                                "5-ring sp2 C", "5-ring HC-N C", "5-ring H2C-N C",
                                "5-ring CH C", "5-ring CH2 C", "5-ring CC4 C",
                                "6-ring aromatic NC=N C", "6-ring aromatic CN C",
                                "6-ring aromatic CF C", "6-ring aromatic with C=O C"]

        bad_hydrogens = ["sp CH H", "sp2 CH2 H", "sp2 CHR H",
                         "sp3 CH H", "sp3 CH2 H", "sp3 CH3 H"]

        bad_nitrogens = ["amine N", "amine NH2 N", "amine NHR N", "amine NRR N",
                         "methylamine N", "dimethylamine N", "trimethylamine N"]

        # second round, reassigns some classifications for GAFF
        for atom in self.atoms.atoms:

            if "ff class done" in atom.tags:
                continue

            if atom.classification in bad_aromatic_carbons:

                if any(c.is_heterocycle for c in atom.cycles):  # cc-cd
                    hets = ["N", "O", "P", "S"]  # etc
                    if any((n.type.real in hets) and
                           ((n.hybridization == "sp3") or
                            n.resonance_sp2) for n in atom.neighbors):
                        atom.type = "cc"
                    else:
                        atom.type = "cd"

                else:
                    atom.type = "ca"

            elif atom.classification == "biphenyl C":  # cp-cq

                if sum(n.classification == "biphenyl C" for
                       n in atom.neighbors) > 1:
                    continue

                while True:
                    if atom.type in ["cp", "cq"]:
                        break
                    atom.type = "cp"
                    for n1 in atom.neighbors:
                        if n1.cycles[0] is not atom.cycles[0]:
                            n1.type = "cp"
                            for n2 in n1.neighbors:
                                if n2.type in ["cp", "cq"]:
                                    continue
                                if n2.classification == "biphenyl C":
                                    n2.type = "cq"
                                    for n3 in n2.neighbors:
                                        if n3.cycles[0] is not n2.cycles[0]:
                                            n3.type = "cq"
                                            for n4 in n3.neighbors:
                                                if n4.type in ["cp", "cq"]:
                                                    continue
                                                if n4.classification == "biphenyl C":
                                                    atom = n4

                # if not any(n.type == "cp" for n in atom.neighbors):
                #     atom.type = "cp"
                # elif (any(n.type == "cp" for n in atom.neighbors)
                #       and not any(any(nn.type == "cp" for nn in n.neighbors) for
                #                   n in atom.neighbors)):
                #     atom.type = "cp"
                # elif not any(n.type == "cq" for n in atom.neighbors):
                #     atom.type = "cq"
                # elif (any(n.type == "cq" for n in atom.neighbors)
                #       and not any(any(nn.type == "cq" for nn in n.neighbors) for
                #                   n in atom.neighbors)):
                #     atom.type = "cq"
                # else:
                #     atom.type = "cp"  # redundant i guess

            elif atom.classification in bad_carbons:

                if atom.hybridization == "sp3":
                    atom.type = "c3"

                elif atom.hybridization == "sp2":

                    if "conjugated" in atom.classification:  # ce-cf

                        if not any(n.type == "c2" for
                                   n in atom.neighbors):
                            continue

                        while True:
                            if atom.type in ["ce", "cf"]:
                                break
                            atom.type = "ce"
                            for n1 in atom.neighbors:
                                if ("conjugated" in n1.classification and
                                        n1.hybridization == "sp2" and
                                        n1.type == "C"):
                                    n1.type = "ce"
                                    for n2 in n1.neighbors:
                                        if ("conjugated" in n2.classification and
                                                n2.hybridization == "sp2" and
                                                n2.type == "C"):
                                            n2.type = "cf"
                                            for n3 in n2.neighbors:
                                                if ("conjugated" in n3.classification and
                                                        n3.hybridization == "sp2" and
                                                        n3.type == "C"):
                                                    n3.type = "cf"
                                                    for n4 in n3.neighbors:
                                                        if ("conjugated" in n4.classification and
                                                                n4.hybridization == "sp2" and
                                                                n4.type == "C"):
                                                            atom = n4

                        # if not any(n.type == "ce" for n in atom.neighbors):
                        #     atom.type = "ce"
                        # elif (any(n.type == "ce" for n in atom.neighbors)
                        #       and not any(any(nn.type == "ce" for nn in n.neighbors) for
                        #                   n in atom.neighbors)):
                        #     atom.type = "ce"
                        # elif not any(n.type == "cf" for n in atom.neighbors):
                        #     atom.type = "cf"
                        # elif (any(n.type == "cf" for n in atom.neighbors)
                        #       and not any(any(nn.type == "cf" for nn in n.neighbors) for
                        #                   n in atom.neighbors)):
                        #     atom.type = "cf"
                        # else:
                        #     atom.type = "ce"  # redundant i guess

                    else:
                        atom.type = "c2"

                elif atom.hybridization == "sp":

                    if "conjugated" in atom.classification:  # cg-ch

                        if not any(n.type == "c1" for
                                   n in atom.neighbors):
                            continue

                        while True:
                            if atom.type in ["cg", "ch"]:
                                break
                            atom.type = "cg"
                            for n1 in atom.neighbors:
                                if ("conjugated" in n1.classification and
                                        n1.hybridization == "sp" and
                                        n1.type == "C"):
                                    n1.type = "cg"
                                    for n2 in n1.neighbors:
                                        if ("conjugated" in n2.classification and
                                                n2.hybridization == "sp" and
                                                n2.type == "C"):
                                            n2.type = "ch"
                                            for n3 in n2.neighbors:
                                                if ("conjugated" in n3.classification and
                                                        n3.hybridization == "sp" and
                                                        n3.type == "C"):
                                                    n3.type = "ch"
                                                    for n4 in n3.neighbors:
                                                        if ("conjugated" in n4.classification and
                                                                n4.hybridization == "sp" and
                                                                n4.type == "C"):
                                                            atom = n4

                        # if not any(n.type == "cg" for n in atom.neighbors):
                        #     atom.type = "cg"
                        # elif (any(n.type == "cg" for n in atom.neighbors)
                        #       and not any(any(nn.type == "cg" for nn in n.neighbors) for
                        #                   n in atom.neighbors)):
                        #     atom.type = "cg"
                        # elif not any(n.type == "ch" for n in atom.neighbors):
                        #     atom.type = "ch"
                        # elif (any(n.type == "ch" for n in atom.neighbors)
                        #       and not any(any(nn.type == "ch" for nn in n.neighbors) for
                        #                   n in atom.neighbors)):
                        #     atom.type = "ch"
                        # else:
                        #     atom.type = "cg"  # redundant i guess

                    else:
                        atom.type = "c1"

            elif atom.classification in bad_hydrogens:

                c = atom.get_neighbors("C")[0]

                # for cases that compete with charmm
                if "aromatic" in c.classification:
                    atom.type = "ha"
                    continue

                electrwd = 0
                for n in c.neighbors:
                    if n.type.real in ["N", "O", "F", "Cl", "Br"]:
                        electrwd += 1

                if "sp3" in atom.classification:
                    if electrwd == 0:
                        atom.type = "hc"
                    elif electrwd == 1:
                        atom.type = "h1"
                    elif electrwd == 2:
                        atom.type = "h2"
                    elif electrwd == 3:
                        atom.type = "h3"

                else:  # sp2 or sp
                    if electrwd == 0:
                        atom.type = "hc"
                    elif electrwd == 1:
                        atom.type = "h4"
                    elif electrwd == 2:
                        atom.type = "h5"

            elif atom.classification in bad_nitrogens:
                if ":S1" in atom.topological_tags:
                    atom.type = "n"  # amide-like e.g. ex. a10
                else:
                    atom.type = "n3"

            elif atom.classification == "imine N":
                if all(n.hybridization == "sp2" for n in atom.neighbors):
                    # conjugated
                    atom.type = "ne"
                    for c in atom.get_neighbors("C"):
                        if len(c.get_neighbors("H")) < 2:  # not terminal
                            # conjugated
                            c.type = "ce"
                else:
                    atom.type = "n2"

        for atom in self.atoms.atoms:

            if "ff class done" in atom.tags:
                continue

            if atom.classification == "5-ring planar 2-bond N":
                for n1 in atom.neighbors:
                    if n1.type == "cd":
                        for n2 in n1.neighbors:
                            if n2.type == "cd":
                                atom.type = "nc"
                if atom.type != "nc":
                    atom.type = "nd"

            elif atom.classification == "6-ring 3-bond N":
                for n1 in atom.neighbors:
                    if n1.type == "ca":
                        atom.type = "nb"
                if atom.type != "nb":
                    atom.type = "na"

            elif atom.classification == "conjugated P":
                for n1 in atom.neighbors:
                    if n1.type == "ce":
                        for n2 in n1.neighbors:
                            if n2.type == "ce":
                                atom.type = "pf"
                if atom.type != "pf":
                    atom.type = "pe"

            elif atom.classification == "aromatic P":
                for n1 in atom.neighbors:
                    if n1.type == "ca":
                        atom.type = "pb"
                    elif n1.type == "cd":
                        for n2 in n1.neighbors:
                            if n2.type == "cd":
                                atom.type = "pc"
                if atom.type not in ["pb", "pc"]:
                    atom.type = "pd"

        # pools
        self.atoms.pool_topological_types(re_compute_types=True)
        self.get_params()

        self.check_water_angles()  # bugs if molecules are periodic!
        # do we use tip3p here or what?

    def get_params(self):
        """Goes through the parameters file, assigning parameter
        values to previously computed topological structures."""
        self.absolute_path = self.par_file
        self.read_file()

        # finds keywords
        atoms_index, bonds_index, angles_index = None, None, None
        dihedrals_index, impropers_index, nonbonded_index = None, None, None
        end_index = None
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
            elif line.startswith("END"):
                end_index = index
            else:
                continue
        assert bonds_index > atoms_index
        assert angles_index > bonds_index
        assert dihedrals_index > angles_index
        assert impropers_index > dihedrals_index
        assert nonbonded_index > impropers_index
        assert end_index > nonbonded_index

        # below, looks for Topological Types that exist in self.atoms
        # and edits their properties based on the parameters file

        # reads ATOMS parameters
        for line in self.content[atoms_index:bonds_index]:

            try:
                atom, mass, *etc = tuple(line.split())
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
                bond = line[0:5]
                k, r0, *etc = tuple(line[5:].split())
            except ValueError:
                continue

            bond = bond.replace("-", ":").replace(" ", "")
            try:
                bond_type = BondType.instances_dict[bond]
            except KeyError:
                try:  # redundant
                    bond_type = BondType.instances_dict[bond]
                except KeyError:
                    continue

            bond_type.k = float(k)
            bond_type.r0 = float(r0)

        # reads ANGLES parameters
        for line in self.content[angles_index:dihedrals_index]:

            try:
                angle = line[0:8]
                k, theta0, *etc = tuple(line[8:].split())
            except ValueError:
                continue

            angle = angle.replace("-", ":").replace(" ", "")
            try:
                angle_type = AngleType.instances_dict[angle]
            except KeyError:
                try:  # redundant
                    angle_type = AngleType.instances_dict[angle]
                except KeyError:
                    continue

            angle_type.k = float(k)
            angle_type.theta0 = float(theta0)

        # reads DIHEDRALS parameters
        for line in self.content[dihedrals_index:impropers_index]:

            try:
                dihedral = line[0:11]
                division, k_term, d, n, *comment = tuple(line[11:].split())
                k = round(float(k_term) / float(division), 5)
                n = int(float(n))
                d = int(float(d))
            except ValueError:
                continue

            dihedral = dihedral.replace("-", ":").replace(" ", "")
            try:
                dihedral_type = DihedralType.instances_dict[dihedral]
            except KeyError:
                continue  # wildcards are checked on a second round
            else:
                try:
                    if not (k in dihedral_type.k and
                            n in dihedral_type.n and
                            d in dihedral_type.d):
                        # to avoid repeated, e.g. h-c-c-c and c-c-c-h
                        dihedral_type.k.append(k)
                        dihedral_type.n.append(n)
                        dihedral_type.d.append(d)
                except (AttributeError, TypeError):
                    dihedral_type.k = [k]
                    dihedral_type.n = [int(float(n))]
                    dihedral_type.d = [int(float(d))]

        # wildcards
        for line in self.content[dihedrals_index:impropers_index]:

            try:
                dihedral = line[0:11]
                division, k_term, d, n, *comment = tuple(line[11:].split())
                k = round(float(k_term) / float(division), 5)
            except ValueError:
                continue

            dihedral = dihedral.replace("-", ":").replace(" ", "")

            dihedral_types = wildcard(DihedralType.instances_dict, dihedral)
            for dihedral_type in dihedral_types:
                dihedral_type.wildcard = True
                try:
                    dihedral_type.k.append(k)
                    dihedral_type.n.append(int(float(n)))
                    dihedral_type.d.append(int(float(d)))
                except AttributeError:
                    dihedral_type.k = [k]
                    dihedral_type.n = [int(float(n))]
                    dihedral_type.d = [int(float(d))]

        # reads IMPROPERS parameters
        for line in self.content[impropers_index:nonbonded_index]:
            # the third atom is the central

            # TODO

            try:
                improper = line[0:11]
                k, phase, v_term, *comment = tuple(line[11:].split())
                assert phase == "180."
            except ValueError:
                continue

            improper = improper.replace("-", ":").replace(" ", "")
            try:
                improper_type = ImproperType.instances_dict[improper]
            except KeyError:
                improper_types = wildcard(ImproperType.instances_dict, improper)
                for improper_type in improper_types:
                    improper_type.k = float(k)
                    improper_type.n = -1  # same as phase = 180.0
            else:
                improper_type.k = float(k)
                improper_type.n = -1  # same as phase = 180.0

        # reads NONBONDED parameters
        r_min_to_sigma = (1 / 2) ** (1 / 6)
        for line in self.content[nonbonded_index:end_index]:

            try:
                atom, r_min_over2, epsilon = tuple(line.split())
            except ValueError:
                continue

            try:
                atom_type = AtomType.instances_dict[atom]
            except KeyError:
                continue

            atom_type.epsilon = round(float(epsilon), 5)
            atom_type.sigma = round(float(r_min_over2) * 2 * r_min_to_sigma, 5)


class Compass(ForceField):
    """"""

    def __init__(self, xyz_file=None):
        super().__init__(xyz_file=xyz_file)
        self.tags.add("compass")
        self.par_file = PACKAGE_DIR + "/compass.prm"
        # TODO if I ever come back to compass, I must use Moltemplate pars

    def compute_topology(self, periodic="", simple=False):

        self.atoms.compute_topology(periodic=periodic, complete=True,
                                    hold_pool_top_types=True,
                                    simple=simple)

        for atom in self.atoms.atoms:
            try:
                # print(atom.classification)
                atom.type = FROM_CLASSIFICATION_TO_COMPASS[atom.classification]
            except KeyError:
                print("WARNING: bad classification for COMPASS: '{}' with [{}]".format(
                    atom.classification, " ".join(list(atom.topological_tags))))

        # pools
        self.atoms.pool_topological_types(re_compute_types=True)
        self.get_params()

    def get_params(self):

        self.absolute_path = self.par_file
        self.read_file()

        # finds keywords
        atoms_index, bonds_index, angles_index = None, None, None
        dihedrals_index, impropers_index, nonbonded_index = None, None, None
        bondbond_index, bondangle_index, angleangle_index = None, None, None
        endbondtorsion_index, middlebondtorsion_index = None, None
        angletorsion_index, angleangletorsion_index = None, None
        for (index, line) in enumerate(self.content):
            if line.startswith("ATOM TYPES"):
                atoms_index = index
            elif line.startswith("NONBONDED"):
                nonbonded_index = index
            elif line.startswith("BONDS"):
                bonds_index = index
            elif line.startswith("ANGLES"):
                angles_index = index
            elif line.startswith("TORSION"):
                dihedrals_index = index
            elif line.startswith("OUT OF PLANE"):
                impropers_index = index
            elif line.startswith("BOND BOND"):
                bondbond_index = index
            elif line.startswith("BOND ANGLE"):
                bondangle_index = index
            elif line.startswith("ANGLE ANGLE"):
                angleangle_index = index
            elif line.startswith("END BOND TORSION"):
                endbondtorsion_index = index
            elif line.startswith("MIDDLE BOND TORSION"):
                middlebondtorsion_index = index
            elif line.startswith("ANGLE TORSION"):
                angletorsion_index = index
            elif line.startswith("ANGLE ANGLE TORSION"):
                angleangletorsion_index = index
            else:
                continue

        # below, looks for Topological Types that exist in self.atoms
        # and edits their properties based on the parameters file

        # reads NONBONDED parameters
        for line in self.content[nonbonded_index:]:

            if line.startswith("END OF SECTION"):
                break

            try:
                atom, r0, e0, *equivalence = tuple(line.split())
            except ValueError:
                continue

            try:
                atom_type = AtomType.instances_dict[atom]
            except KeyError:
                continue

            atom_type.epsilon = float(e0)
            atom_type.sigma = float(r0)

        # reads BONDS parameters
        for line in self.content[bonds_index:]:

            if line.startswith("END OF SECTION"):
                break

            try:
                atom1, atom2, b0, k2, k3, k4, increment = tuple(line.split())
            except ValueError:
                continue

            #if atom1 > atom2:
            #    atom1, atom2 = atom2, atom1
            bond = atom1 + ":" + atom2
            try:
                bond_type = BondType.instances_dict[bond]
            except KeyError:
                try:  # redundant
                    bond_type = BondType.instances_dict[bond]
                except KeyError:
                    continue

            bond_type.r0 = float(b0)
            bond_type.k = float(k2)
            bond_type.k3 = float(k3)
            bond_type.k4 = float(k4)
            bond_type.increment = float(increment)
            bond_type.first_atom = atom1

        # reads ANGLES parameters
        for line in self.content[angles_index:]:

            if line.startswith("END OF SECTION"):
                break

            try:
                atom1, atom2, atom3, q0, k2, k3, k4 = tuple(line.split())
            except ValueError:
                continue

            if atom1 > atom3:
                atom1, atom3 = atom3, atom1

            angle = atom1 + ":" + atom2 + ":" + atom3
            try:
                angle_type = AngleType.instances_dict[angle]
            except KeyError:
                try:  # redundant
                    angle_type = AngleType.instances_dict[angle]
                except KeyError:
                    continue

            angle_type.theta0 = float(q0)
            angle_type.k = float(k2)
            angle_type.k3 = float(k3)
            angle_type.k4 = float(k4)

            bond1_type = BondType.instances_dict[atom1 + ":" + atom2]
            angle_type.r1 = bond1_type.r0
            bond2_type = BondType.instances_dict[atom2 + ":" + atom3]
            angle_type.r2 = bond2_type.r0

        # reads DIHEDRALS parameters
        for line in self.content[dihedrals_index:]:

            if line.startswith("END OF SECTION"):
                break

            try:
                atom1, atom2, atom3, atom4, k1, k2, k3 = tuple(line.split())
            except ValueError:
                continue

            if atom1 > atom4:
                atom1, atom2, atom3, atom4 = atom4, atom3, atom2, atom1
            elif (atom1 == atom4) and (atom2 > atom3):
                atom1, atom2, atom3, atom4 = atom4, atom3, atom2, atom1

            dihedral = atom1 + ":" + atom2 + ":" + atom3 + ":" + atom4
            try:
                dihedral_type = DihedralType.instances_dict[dihedral]
            except KeyError:
                try:  # redundant
                    dihedral_type = DihedralType.instances_dict[dihedral]
                except KeyError:
                    continue

            dihedral_type.k = float(k1)
            dihedral_type.k2 = float(k2)
            dihedral_type.k3 = float(k3)

            bond1_type = BondType.instances_dict[atom1 + ":" + atom2]
            dihedral_type.r1 = bond1_type.r0
            bond2_type = BondType.instances_dict[atom2 + ":" + atom3]
            dihedral_type.r2 = bond2_type.r0
            bond3_type = BondType.instances_dict[atom3 + ":" + atom4]
            dihedral_type.r3 = bond3_type.r0
            angle1_type = AngleType.instances_dict[atom1 + ":" + atom2 + ":" + atom3]
            dihedral_type.theta1 = angle1_type.theta0
            angle2_type = AngleType.instances_dict[atom2 + ":" + atom3 + ":" + atom4]
            dihedral_type.theta2 = angle2_type.theta0

        # reads IMPROPERS parameters
        for line in self.content[impropers_index:]:

            if line.startswith("END OF SECTION"):
                break

            try:
                atom1, atom2, atom3, atom4, k2 = tuple(line.split())
            except ValueError:
                continue

            improper = atom2 + ":" + atom1 + ":" + atom3 + ":" + atom4
            # the second atom is the central in COMPASS!
            try:
                improper_type = ImproperType.instances_dict[improper]
            except KeyError:
                continue

            improper_type.k = float(k2)

        # reads BOND BOND parameters
        for line in self.content[bondbond_index:]:

            if line.startswith("END OF SECTION"):
                break

            try:
                i, j, k, ij_jk = tuple(line.split())
            except ValueError:
                continue

            angle = i + ":" + j + ":" + k
            try:
                angle_type = AngleType.instances_dict[angle]
            except KeyError:
                continue

            angle_type.ij_jk = float(ij_jk)

        # reads BOND ANGLE parameters
        for line in self.content[bondangle_index:]:

            if line.startswith("END OF SECTION"):
                break

            try:
                i, j, k, ij_ijk, jk_ijk = tuple(line.split())
            except ValueError:
                continue

            angle = i + ":" + j + ":" + k
            try:
                angle_type = AngleType.instances_dict[angle]
            except KeyError:
                continue

            angle_type.ij_ijk = float(ij_ijk)
            angle_type.jk_ijk = float(jk_ijk)

        # reads ANGLE ANGLE parameters
        # (these go as Improper pars)
        for line in self.content[angleangle_index:]:
            # TODO this is wrong
            #  and weird and I need to look this up if I ever need AngleAngle

            if line.startswith("END OF SECTION"):
                break

            try:
                i, j, k, l, ijk_jkl = tuple(line.split())
                # there's actually 3 energy terms; see LAMMPS doc for this
            except ValueError:
                continue

            improper = j + ":" + i + ":" + k + ":" + l
            # the second atom is the central in COMPASS!
            try:
                improper_type = ImproperType.instances_dict[improper]
            except KeyError:
                continue

            # improper_type.ijk_jkl = float(ijk_jkl)

        # reads END BOND TORSION parameters
        for line in self.content[endbondtorsion_index:]:

            if line.startswith("END OF SECTION"):
                break

            try:
                i, j, k, l, b1, b2, b3, c1, c2, c3 = tuple(line.split())
            except ValueError:
                continue

            dihedral = i + ":" + j + ":" + k + ":" + l
            try:
                dihedral_type = DihedralType.instances_dict[dihedral]
            except KeyError:
                continue

            dihedral_type.b1 = float(b1)
            dihedral_type.b2 = float(b2)
            dihedral_type.b3 = float(b3)
            dihedral_type.c1 = float(c1)
            dihedral_type.c2 = float(c2)
            dihedral_type.c3 = float(c3)

        # reads MIDDLE BOND TORSION parameters
        for line in self.content[middlebondtorsion_index:angletorsion_index]:

            if line.startswith("END OF SECTION"):
                break

            try:
                i, j, k, l, a1, a2, a3 = tuple(line.split())
            except ValueError:
                continue

            dihedral = i + ":" + j + ":" + k + ":" + l
            try:
                dihedral_type = DihedralType.instances_dict[dihedral]
            except KeyError:
                continue

            dihedral_type.a1 = float(a1)
            dihedral_type.a2 = float(a2)
            dihedral_type.a3 = float(a3)

        # reads ANGLE TORSION parameters
        for line in self.content[angletorsion_index:angleangletorsion_index]:

            if line.startswith("END OF SECTION"):
                break

            try:
                i, j, k, l, d1, d2, d3, e1, e2, e3 = tuple(line.split())
            except ValueError:
                continue

            dihedral = i + ":" + j + ":" + k + ":" + l
            try:
                dihedral_type = DihedralType.instances_dict[dihedral]
            except KeyError:
                continue

            dihedral_type.d1 = float(d1)
            dihedral_type.d2 = float(d2)
            dihedral_type.d3 = float(d3)
            dihedral_type.e1 = float(e1)
            dihedral_type.e2 = float(e2)
            dihedral_type.e3 = float(e3)

        # reads ANGLE ANGLE TORSION parameters
        for line in self.content[angleangletorsion_index:]:

            if line.startswith("END OF SECTION"):
                break

            try:
                i, j, k, l, ij_ijkl_jkl = tuple(line.split())
            except ValueError:
                continue

            dihedral = i + ":" + j + ":" + k + ":" + l
            try:
                dihedral_type = DihedralType.instances_dict[dihedral]
            except KeyError:
                continue

            dihedral_type.ij_ijkl_jkl = float(ij_ijkl_jkl)

        # TODO:
        # BondBond13 ???

    def order_in_bonds(self):
        for bond in self.atoms.bonds:
            atom1, atom2 = tuple(bond.atoms)
            if atom1.type != bond.type.first_atom:
                bond.atoms = [atom2, atom1]
        for bond_type in self.atoms.bond_types:
            bond_type.class2_order()

    def increment_charges(self):

        for atom in self.atoms:
            atom.charge = 0.0000

        for bond in self.atoms.bonds:
            bond_type = bond.type
            increment = round(bond_type.increment, 4)
            # print(bond_type, increment, bond_type.first_atom)
            for atom in bond:
                if str(atom.type) == str(bond_type.first_atom):
                    # print("Adding to {}".format(atom.type))
                    atom.charge += increment
                else:
                    # print("Taking from {}".format(atom.type))
                    atom.charge -= increment

        for atom in self.atoms:
            atom.charge = round(atom.charge, 4)

        # checks neutrality
        total_charge = 0.0
        for atom in self.atoms:
            total_charge += atom.charge
        print("System total charge = {}".format(total_charge))

    def order_in_angles(self):
        """Makes sure every angle has atoms ordered in alphabetical order
        (using their COMPASS types). This ensures that r1 and r2 correspond
        to the right bonds, in order, which is needed for class2 parameters."""
        for angle in self.atoms.angles:
            atom1, atom2, atom3 = tuple(angle.atoms)
            if str(atom1.type) > str(atom3.type):
                angle.atoms = [atom3, atom2, atom1]

    def order_in_dihedrals(self):
        """Same as order_in_angles, but for dihedrals."""
        for dihedral in self.atoms.dihedrals:
            atom1, atom2, atom3, atom4 = tuple(dihedral.atoms)
            if str(atom1.type) > str(atom4.type):
                dihedral.atoms = [atom4, atom3, atom2, atom1]
            elif (str(atom1.type) == str(atom4.type) and
                  str(atom2.type) > str(atom3.type)):
                dihedral.atoms = [atom4, atom3, atom2, atom1]

    def write_par(self, output_path):
        """Check docstring of ForceField.write_par()"""

        with open(output_path, "w") as F:

            F.write("Atom Types\n\n")
            F.write("# type\t'charge'\tepsilon\tsigma\n")
            for atom_type in self.atoms.atom_types:
                F.write("\t".join([str(atom_type), str(0.0000),
                                   str(atom_type.epsilon), str(atom_type.sigma)]) + "\n")

            F.write("\nBond Types\n\n")
            F.write("# bond_type\tincrement\tr0\tk2\tk3\tk4\n")
            for bond_type in self.atoms.bond_types:
                F.write("\t".join([str(bond_type), str(bond_type.increment),
                                   str(bond_type.r0), str(bond_type.k),
                                   str(bond_type.k3), str(bond_type.k4)]) + "\n")

            F.write("\nAngle Types\n\n")
            F.write("# angle_type\ttheta0\tk2\tk3\tk4\n")
            for angle_type in self.atoms.angle_types:
                F.write("\t".join([str(angle_type),
                                   str(angle_type.theta0), str(angle_type.k),
                                   str(angle_type.k3), str(angle_type.k4)]) + "\n")

            F.write("\nDihedral Types\n\n")
            F.write("# dihedral_type\tk1\tphi1\tk2\tphi2\tk3\tphi3\n")
            for dihedral_type in self.atoms.dihedral_types:
                F.write("\t".join([str(dihedral_type),
                                   str(dihedral_type.k), str(dihedral_type.phi1),
                                   str(dihedral_type.k2), str(dihedral_type.phi2),
                                   str(dihedral_type.k3), str(dihedral_type.phi3)]) + "\n")

            F.write("\nImproper Types\n\n")
            F.write("# improper_type\tk\tx0=0\n")
            for improper_type in self.atoms.improper_types:
                F.write("\t".join([improper_type.compass_order,
                                   str(improper_type.k), str(improper_type.x0)]) + "\n")

            F.write("\nBondBond Coeffs\n\n")
            F.write("# angle_type\tm\tr1\tr2\n")
            for angle_type in self.atoms.angle_types:
                F.write("\t".join([str(angle_type), str(angle_type.ij_jk), str(angle_type.r1),
                                   str(angle_type.r2)]) + "\n")

            F.write("\nBondAngle Coeffs\n\n")
            F.write("# angle_type\tn1\tn2\tr1\tr2\n")
            for angle_type in self.atoms.angle_types:
                F.write("\t".join([str(angle_type), str(angle_type.ij_ijk), str(angle_type.jk_ijk),
                                   str(angle_type.r1), str(angle_type.r2)]) + "\n")

            F.write("\nMiddleBondTorsion Coeffs\n\n")
            F.write("# dihedral_type\ta1\ta2\ta3\tr2\n")
            for dihedral_type in self.atoms.dihedral_types:
                F.write("\t".join([str(dihedral_type),
                                   str(dihedral_type.a1), str(dihedral_type.a2),
                                   str(dihedral_type.a3), str(dihedral_type.r2)]) + "\n")

            F.write("\nEndBondTorsion Coeffs\n\n")
            F.write("# dihedral_type\tb1\tb2\tb3\tc1\tc2\tc3\tr1\tr3\n")
            for dihedral_type in self.atoms.dihedral_types:
                F.write("\t".join([str(dihedral_type),
                                   str(dihedral_type.b1), str(dihedral_type.b2),
                                   str(dihedral_type.b3), str(dihedral_type.c1),
                                   str(dihedral_type.c2), str(dihedral_type.c3),
                                   str(dihedral_type.r1), str(dihedral_type.r3)]) + "\n")

            F.write("\nAngleTorsion Coeffs\n\n")
            F.write("# dihedral_type\td1\td2\td3\te1\te2\te3\ttheta1\ttheta2\n")
            for dihedral_type in self.atoms.dihedral_types:
                F.write("\t".join([str(dihedral_type),
                                   str(dihedral_type.d1), str(dihedral_type.d2),
                                   str(dihedral_type.d3), str(dihedral_type.e1),
                                   str(dihedral_type.e2), str(dihedral_type.e3),
                                   str(dihedral_type.theta1), str(dihedral_type.theta2)]) + "\n")

            F.write("\nAngleAngleTorsion Coeffs\n\n")
            F.write("# dihedral_type\tij_ijkl_jkl\ttheta1\ttheta2\n")
            for dihedral_type in self.atoms.dihedral_types:
                F.write("\t".join([str(dihedral_type),
                                   str(dihedral_type.ij_ijkl_jkl), str(dihedral_type.theta1),
                                   str(dihedral_type.theta2)]) + "\n")

            F.write("\nBondBond13 Coeffs\n\n")
            F.write("# dihedral_type\tn\tr1\tr3\n")
            for dihedral_type in self.atoms.dihedral_types:
                F.write("\t".join([str(dihedral_type),
                                   str(dihedral_type.n), str(dihedral_type.r1),
                                   str(dihedral_type.r3)]) + "\n")

            F.write("\nAngleAngle Coeffs\n\n")
            F.write("# improper_type\tm1\tm2\tm3\ttheta1\ttheta2\ttheta3\n")
            for improper_type in self.atoms.improper_types:
                F.write("\t".join([str(improper_type),
                                   str(improper_type.m1), str(improper_type.m2),
                                   str(improper_type.m3), str(improper_type.theta1),
                                   str(improper_type.theta2), str(improper_type.theta3)]) + "\n")
