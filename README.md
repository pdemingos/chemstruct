# ChemStruct - Chemical Analysis Package

A package for structural analysis of atomic files (e.g. xyz files)
as well as for assistance on molecular dynamics simulations. 

Its main objective is to take an xyz file with nothing but atomic species
and positions (and a 3x3 lattice matrix if the system is periodic) 
and output a LAMMPS Data File ready for simulation, with the atoms 
classified according to a force field and with the force field parameters
already set for the existing atoms, bonds, angles and dihedrals. 
The current version has the CHARMM General Force Field implemented in that way.

Ready-to-run examples can be found in the directory with that name. 

