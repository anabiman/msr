/*
 * @Author: Andrew Abi Mansour
 * @Created: March Jan 14, 2013
 *
 * MSR is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * MSR is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with DMS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 */

#include "main.h"
static char help[] = "Atomic Sparse Reconstruction\n";

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char** argv) {

  PetscInitialize(&argc, &argv, NULL, help);

  // System-dependent params
  const int dims(3);
  PetscInt numCoords, numCG, numCons;
  PetscReal tol;

  // I/O file names
  char coordFname[PETSC_MAX_PATH_LEN],
       indicesFname[PETSC_MAX_PATH_LEN],
       eqLengthFname[PETSC_MAX_PATH_LEN],
       oFname[PETSC_MAX_PATH_LEN],
       refFname[PETSC_MAX_PATH_LEN],
       phiFname[PETSC_MAX_PATH_LEN];

  PetscBool flg = readInput(coordFname, phiFname, indicesFname, eqLengthFname, refFname, oFname, numCoords, numCons, numCG);
  PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "--tol", &tol, &flg);

  // mpi params
  PetscInt numProcs, rank;
  MPI_Comm_size(PETSC_COMM_WORLD, &numProcs);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  Vec CoordsX, CoordsY, CoordsZ, PhiX, PhiY, PhiZ;
  Mat Ref;
  Vec EqLength;

  PetscScalar X, Y, Z;
  PetscErrorCode ierr;

  PetscInt *indicesOne,
	   *indicesTwo;

  VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, numCoords, &CoordsX);
  VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, numCoords, &CoordsY);
  VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, numCoords, &CoordsZ);

  VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, numCG, &PhiX);
  VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, numCG, &PhiY);
  VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, numCG, &PhiZ);

  if(rank >= 0){  // read input on all processors

    if (flg) {
          	readMatrixASCII(Ref, refFname, numCoords, numCG, "%le");
          	readVectorASCII(EqLength, eqLengthFname, numCons,"%le");

		PetscInt istart, iend;
  		VecGetOwnershipRange(EqLength,&istart,&iend);

		indicesOne = new PetscInt[iend-istart],
	 	indicesTwo = new PetscInt[iend-istart];

        	std::ifstream fp;
                fp.open(indicesFname, std::ios::in);

		PetscInt tmp;

        	for(auto i = 0; i < numCons; i++)
			if(i >= istart && i < iend)
        			fp >> indicesOne[i - istart] >> indicesTwo[i - istart]; 
			else
				if(i < iend - 1)
					fp >> tmp >> tmp;

        	fp.close();

        	fp.open(coordFname, std::ios::in);
		fp.clear();
		fp.seekg(0, std::ios::beg);

        	for(auto i = 0; i < numCoords; i++) {

            		fp >> X >> Y >> Z;

        		VecSetValues(CoordsX, 1, &i, &X, INSERT_VALUES);
        		VecSetValues(CoordsY, 1, &i, &Y, INSERT_VALUES);
        		VecSetValues(CoordsZ, 1, &i, &Z, INSERT_VALUES);
        	}

	         fp.close();

	          fp.open(phiFname, std::ios::in);

            for(auto i = 0; i < numCG; i++) {

                fp >> X >> Y >> Z;

                VecSetValues(PhiX, 1, &i, &X, INSERT_VALUES);
                VecSetValues(PhiY, 1, &i, &Y, INSERT_VALUES);
                VecSetValues(PhiZ, 1, &i, &Z, INSERT_VALUES);
        }

    }
  }

    std::vector<Vec> Coords(dims);
    Coords[0] = CoordsX;
    Coords[1] = CoordsY;
    Coords[2] = CoordsZ;

    std::vector<Vec> Phi(dims);
    Phi[0] = PhiX;
    Phi[1] = PhiY;
    Phi[2] = PhiZ;

    for(auto dim = 0; dim < dims; dim++) {

            ierr = VecAssemblyBegin(Coords[dim]); CHKERRQ(ierr);
	    ierr = VecAssemblyBegin(Phi[dim]); CHKERRQ(ierr);

            ierr = VecAssemblyEnd(Coords[dim]); CHKERRQ(ierr);
	    ierr = VecAssemblyEnd(Phi[dim]); CHKERRQ(ierr);
    }

      /**************************************/
     /****** Begin Computation ~ Phew! *****/
    /**************************************/

    NewtonIter(Coords, Phi, indicesOne, indicesTwo, EqLength, Ref, oFname, tol);
    writeVector(Coords, oFname);

      /**************************************************/
     /****** End Computation ~ Garbage collection *****/
    /*************************************************/

    cleanUp(Ref, EqLength, Coords, Phi);

    PetscFinalize();
    return 0;
}

