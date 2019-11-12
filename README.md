Microstate Sparse Reconstruction
================================
Microstate sparse reconstruction (MSR) is a backmapping algorithm that recovers all-atom structures from coarse-grained/reduced variables. MSR uses topological information (connectivity) of a given macromolecule and methods from optimization and inverse theory to rapidly achieve this recovery.

The schematic below shows how MSR can be used to reconstruct a physically meaningful structure of a strongly perturbed alanine dipeptide. 
<p style="text-align:center;"><img src="data/imgs/recovery.png"></p>

AUTHOR
======
Andrew Abi-Mansour
Department of Chemistry, Indiana University, Bloomington

**Please consider citing the following paper if your find this code useful in your research:**

[![DOI for Citing MSR](https://img.shields.io/badge/DOI-10.1021%2Facs.jctc.5b00056-blue.svg)](https://doi.org/10.1021/acs.jctc.6b00348)

PREREQUISITES
=============
PETSc - https://www.mcs.anl.gov/petsc

MDAnalysis - https://code.google.com/p/mdanalysis

Numpy - http://www.numpy.org

OpenMPI - http://www.open-mpi.org or MPICH  - https://www.mpich.org

PYTHON COMMAND-LINE INTERFACE
=============================
* --top: input reference structure (topology) file (pdb, gro, ...) 
* --cgOP: coarse-graining matrix file (table) or python code (module)
* --fgOP: fine-graining (backmapping) matrix file (table) or python code (module)
* [--out]: output filename of recovered all-atom structure 
* [--tol]: tolerance set for the atomic displacement, below which convergence is assumed to be achieved, defaults to 0.1A
* [--maxiter]: max number of iterations the solver performs before giving up, defaults to 100

C++ COMMAND-LINE INTERFACE
==========================

* --ref: path to input reference all-atom coordinates
* --indices: path to input bond/angle indices
* --lengths: path to input bond/angle lengths
* --cg: path to input CG coordinates
* --cgOP: path to the coarse-graining operator (matrix) file
* --fgOP: path to inverse (backmapping) operator (matrix) file
* [--out]: output filename of recovered all-atom positions
* [--tol]: tolerance of the iterative solver

WORKING EXAMPLE
===============
A sample script (swm_example.py) is provided in the examples directory. This script has the following signature:

`python swm_example.py top tol`

* top: topology file such as a pdb or gro file
* tol: tolerance set for the atomic displacement, below which convergence is assumed to be achieved

This script uses the space-warping method to coarse-grained an all-atom system perturbed by noise. MSR recovers a reasonable microstate in 20 iterations.
To run this script, invoke:

`swm_example.py ../data/systems/dialanine.pdb 1e-3`

FEEDBACK
========
Email me: andrew [at] gmail [dot] com
