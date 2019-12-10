%module msr
 %{
 /* Includes the header in the wrapper code */
 #define SWIG_FILE_WITH_INIT
 #include "../src/python_api.h"
 %}

%include "numpy.i"
%init %{
import_array();
%}

%apply (double *ARGOUT_ARRAY1, int DIM1) {(double *COORDS_OUT, int NATOMS_BY_3)};
%apply (double *IN_ARRAY2, int DIM1, int DIM2) {(double *COORDS_IN, int NATOMS2, int DIMC2)};
%apply (double *IN_ARRAY2, int DIM1, int DIM2) {(double *CG_IN, int NCG1, int DIMCG)};
%apply (int *IN_ARRAY2, int DIM1, int DIM2) {(int *INDICES, int NCONS1, int DIMCN)};
%apply (double *IN_ARRAY1, int DIM1) {(double *LENGTHS, int NCONS2)};
%apply (double *IN_ARRAY2, int DIM1, int DIM2) {(double *CGOP, int NCG, int NATOMS)};


// NewtonIter(std::vector<Vec>& Coords, std::vector<Vec>& Phi, PetscInt* indicesOne, PetscInt* indicesTwo, Vec& EqLength, Mat& Ref)

%feature("autodoc", "Computes an all-arom state: reconstruct(coords, CG, indices, lengths, CGOP) -> coords");
%include "../src/python_api.h"

