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

#include "Main.h"
static char help[] = "Atomic Sparse Reconstruction\n";
using namespace std;

#undef __FUNCT__
#define __FUNCT__ "ReadVectorASCII"
PetscErrorCode ReadVectorASCII(Vec& Vector, char* fname, const PetscInt nrows, const char* type) {
  /* Read in vector from ascii files */
  PetscFunctionBegin;
  PetscErrorCode ierr;

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  ierr = VecCreate(PETSC_COMM_WORLD,&Vector); CHKERRQ(ierr);
  ierr = VecSetSizes(Vector,PETSC_DECIDE,nrows); CHKERRQ(ierr);
  ierr = VecSetType(Vector,VECMPI); CHKERRQ(ierr);

  if(!rank) {
  	FILE *fp;
  	ierr = PetscFOpen(PETSC_COMM_SELF,fname,"r",&fp); CHKERRQ(ierr);

  	for (int i = 0; i < nrows; i++) {
      		PetscScalar vals;
      		fscanf(fp,type,&vals);
      		fscanf(fp,"\n");
     		 ierr = VecSetValues(Vector,1,&i,&vals,INSERT_VALUES);CHKERRQ(ierr);
    	}

  fclose(fp);

  }

  VecAssemblyBegin(Vector);
  VecAssemblyEnd(Vector);
  fflush(stdout);

  PetscFunctionReturn(ierr);
}

#undef __FUNCT__
#define __FUNCT__ "ReadMatrixASCII"
PetscErrorCode ReadMatrixASCII(Mat& Matrix, char* fname, const PetscInt nrows, const PetscInt ncols, const char* type) {
  /* Read in matrix from ascii files */
  PetscFunctionBegin;
  PetscErrorCode ierr;

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  ierr = MatCreateDense(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, nrows, ncols, NULL, &Matrix);

  try {

     if(!rank) {
      	FILE *fp;
      	ierr = PetscFOpen(PETSC_COMM_SELF, fname, "r", &fp); CHKERRQ(ierr);

      	for (int i = 0; i < nrows; i++) {
        	PetscScalar vals[ncols];
          	PetscInt cols[ncols];

          	for (int j = 0; j < ncols; j++) {
              		fscanf(fp,type,&vals[j]);
              		cols[j] = j;
          	}

          	fscanf(fp,"\n");
          	ierr = MatSetValues(Matrix,1,&i,ncols,cols,vals,INSERT_VALUES); CHKERRQ(ierr);
        }

      fclose(fp);

      }

      MatAssemblyBegin(Matrix,MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(Matrix,MAT_FINAL_ASSEMBLY);
      //fflush(stdout);

    }
  catch(int e) {
    PetscFinalize();
    std::cout << "Failure in function " << __FUNCT__ << " with error " << e << std::endl;
  }

  PetscFunctionReturn(ierr);
}

#undef __FUNCT__
#define __FUNCT__ "writeVector"
PetscErrorCode writeVector(std::vector<Vec>& Coords, const char* fname) {

  PetscFunctionBegin;
  PetscErrorCode ierr;
  PetscViewer viewer;

  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII); CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_APPEND); CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, fname); CHKERRQ(ierr);

  for(auto it = Coords.begin(); it != Coords.end(); it++)
      ierr = VecView(*it, viewer);

  CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}

#undef __FUNCT__
#define __FUNCT__ "readInput"
PetscBool readInput(char* coordFname, char* phiFname, char* indicesFname, char* eqLengthFname, char* refFname, char* refTrans,
		    char* oFname, PetscInt& numCoords, PetscInt& numCons, PetscInt& numCG)
{

    PetscFunctionBegin;
    PetscBool flg;

    PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "-c", coordFname, sizeof(coordFname), &flg);
    PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "-cg", phiFname, sizeof(phiFname), &flg);
    PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "-i", indicesFname, sizeof(indicesFname), &flg);
    PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "-l", eqLengthFname, sizeof(eqLengthFname), &flg);
    PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "-ref", refFname, sizeof(refFname), &flg);
    PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "-refTrans", refTrans, sizeof(refTrans), &flg);
    PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "-o", oFname, sizeof(oFname), &flg);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-nc", &numCoords, &flg);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-ns", &numCons, &flg);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-ncg", &numCG, &flg);

    PetscFunctionReturn(flg);
}

#undef __FUNCT__
#define __FUNCT__ "cleanUp"
PetscErrorCode cleanUp(Mat& Ref, Mat& RefTrans, Vec& EqLength, std::vector<Vec>& Coords, std::vector<Vec>& Phi ) {

    PetscFunctionBegin;
    PetscErrorCode ierr;

    ierr = MatDestroy(&Ref); CHKERRQ(ierr);
    ierr = MatDestroy(&RefTrans); CHKERRQ(ierr);
    ierr = VecDestroy(&EqLength); CHKERRQ(ierr);

    for(auto it = Coords.begin(); it != Coords.end(); it++) {
      ierr = VecDestroy(&(*it));
      CHKERRQ(ierr);
    }

    for(auto it = Phi.begin(); it != Phi.end(); it++) {
      ierr = VecDestroy(&(*it));
      CHKERRQ(ierr);
    }

    PetscFunctionReturn(ierr);
}

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
       phiFname[PETSC_MAX_PATH_LEN],
       refTrans[PETSC_MAX_PATH_LEN];

  PetscBool flg = readInput(coordFname, phiFname, indicesFname, eqLengthFname, refFname, refTrans, oFname, numCoords, numCons, numCG);
  PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-tol", &tol, &flg);

  // mpi params
  PetscInt numProcs, rank;
  MPI_Comm_size(PETSC_COMM_WORLD, &numProcs);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  Vec CoordsX, CoordsY, CoordsZ, PhiX, PhiY, PhiZ;
  Mat Ref, RefTrans;
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
          	ReadMatrixASCII(Ref, refFname, numCoords, numCG, "%le");
		ReadMatrixASCII(RefTrans, refTrans, numCG, numCG, "%le");
          	ReadVectorASCII(EqLength, eqLengthFname, numCons,"%le");

		PetscInt istart, iend;
  		VecGetOwnershipRange(EqLength,&istart,&iend);

		indicesOne = new PetscInt[iend-istart],
	 	indicesTwo = new PetscInt[iend-istart];

        	ifstream fp;
                fp.open(indicesFname, ios::in);

		PetscInt tmp;
		
        	for(auto i = 0; i < numCons; i++)
			if(i >= istart && i < iend)
        			fp >> indicesOne[i - istart] >> indicesTwo[i - istart]; 
			else
				if(i < iend - 1)
					fp >> tmp >> tmp;

        	fp.close();

        	fp.open(coordFname, ios::in);
		fp.clear();
		fp.seekg(0, ios::beg);

        	for(auto i = 0; i < numCoords; i++) {

            		fp >> X >> Y >> Z;

        		VecSetValues(CoordsX, 1, &i, &X, INSERT_VALUES);
        		VecSetValues(CoordsY, 1, &i, &Y, INSERT_VALUES);
        		VecSetValues(CoordsZ, 1, &i, &Z, INSERT_VALUES);
        	}

	         fp.close();

	          fp.open(phiFname, ios::in);

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

   //VecView(Coords[0], PETSC_VIEWER_STDOUT_SELF);

      /**************************************/
     /****** Begin Computation ~ Phew! *****/
    /**************************************/

    NewtonIter(Coords, Phi, indicesOne, indicesTwo, EqLength, Ref, RefTrans, oFname, tol);
    writeVector(Coords, oFname);

      /**************************************************/
     /****** End Computation ~ Garbage collection *****/
    /*************************************************/

    cleanUp(Ref, RefTrans, EqLength, Coords, Phi);

    PetscFinalize();
    return 0;
}
