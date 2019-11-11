MSR
===
Microstate sparse reconstruction (MSR) is an optimization algorithm that solves the inverse problem: recovering an all-atom state consistent with a coarse-grained description.

**Please consider citing the following paper if your find MSR useful in your research:**

[![DOI for Citing MSR](https://img.shields.io/badge/DOI-10.1021%2Facs.jctc.5b00056-blue.svg)](https://doi.org/10.1021/acs.jctc.6b00348)


AUTHOR
======
Andrew Abi-Mansour

INSTITUTION
===========
Department of Chemistry, Indiana University, Bloomington

PREREQUISITES
=============
ProtoMD - https://github.com/CTCNano/proto_md

PETSc - https://www.mcs.anl.gov/petsc

MDAnalysis - https://code.google.com/p/mdanalysis

Numpy - http://www.numpy.org

OpenMPI - http://www.open-mpi.org or MPICH  - https://www.mpich.org

USAGE
=====
`mpirun MSR.a -nc {int} -ns {int} -c {filename} -i {filename} -l {filename} -cg {filename} -ncg {int} -ref {filename} -inv {filename} -o {filename} -tol {float}`

* nc: number of atoms
* ns: number of constraints
* c: path to input reference all-atom coordinates
* i: path to input bond/angle indices
* l: path to input bond/angle lengths
* cg: path to input CG coordinates
* ref: path to reference basis functions (matrix) file
* inv: path to inverse operator (matrix) file
* o: output filename of recovered all-atom positions
* tol: tolerance of the iterative solver

WORKING EXAMPLE
===============
A sample script (swm_example.py) is provided in the examples directory. This script has the following signature:

`python swm_example.py top tol`

* top: topology file such as a pdb or gro file
* tol: tolerance set for the atomic displacement, below which convergence is assumed to be achieved

This script uses the space-warping method to coarse-grained an all-atom system perturbed by noise. MSR recovers a reasonable microstate in 20 iterations.
To run this script, invoke:

`swm_example.py ../data/systems/dialanine.pdb 1e-3`
