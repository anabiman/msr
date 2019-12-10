#ifndef PYTHON_API
#define PYTHON_API

// double (std::vector<Vec>&, std::vector<Vec>&, PetscInt*, PetscInt*, Vec&, Mat&, Mat&, char*, PetscReal);
//PetscErrorCode NewIter(std::vector<Vec>&, std::vector<Vec>&, PetscInt*, PetscInt*, Vec&, Mat&, Mat&, char*, PetscReal);

void reconstruct(double *COORDS_IN, int NATOMS2, int DIMC2, double *CG_IN, int NCG1, int DIMCG, int* INDICES, 
	int NCONS1, int DIMCN, double* LENGTHS, int NCONS2, double* CGOP, int NCG, int NATOMS, double *COORDS_OUT, int NATOMS_BY_3) {

	for(int i = 0; i < NATOMS_BY_3; i++)
		COORDS_OUT[i] = i;

}

#endif
