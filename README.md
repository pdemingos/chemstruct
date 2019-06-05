ChemStruct: Chemical Structure Analysis Package
================================================

## About

ChemStruct is a Python package for structural analysis of atomic files (like xyz files) as well as for assistance on molecular dynamics simulations. 

Its main objective is to take an xyz file with nothing but atomic species and positions (and a 3x3 lattice matrix if the system is periodic) and output a LAMMPS Data File ready for simulation, with the atoms classified according to a force field and with the force field parameters already set for the existing atoms, bonds, angles and dihedrals. 

Ready-to-run examples can be found in the `examples` directory.

## Usage

Install ChemStruct using pip:

```bash
pip3 install chemstruct
```

Then run Python, import the package and use the `write_lmpdat()` function to convert your xyz file into a LAMMPS Data File. 

```bash
python3
```

```bash
from chemstruct.quick import write_lmpdat
write_lmpdat("yourdir/yourfile.xyz", "outdir/outfile.lmp")
```
