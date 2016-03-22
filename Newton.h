#ifndef NEWTON_H
#define NEWTON_H

#include <iostream>
#include <fstream>
#include <string>
#include <petscksp.h>
#include <petscmat.h>
#include <vector>
#include <algorithm>

PetscErrorCode NewtonIter(std::vector<Vec>&, std::vector<Vec>&, PetscInt*, PetscInt*, Vec&, Mat&, Mat&, char* Fname, PetscReal);
PetscErrorCode writeVector(std::vector<Vec>& Coords, const char* fname);

#endif
