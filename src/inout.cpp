/*
 * @Author: Andrew Abi Mansour
 * @Created: Nov 12, 2019
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

#undef __FUNCT__
#define __FUNCT__ "countLines"
int countLines(std::string filename) {
    int number_of_lines = 0;
    std::string line;
    std::ifstream myfile(filename);

    while (std::getline(myfile, line))
        ++number_of_lines;

    return number_of_lines;
}

#undef __FUNCT__
#define __FUNCT__ "readVectorASCII"
PetscErrorCode readVectorASCII(Vec& Vector, char* fname, const PetscInt nrows, const char* type) {
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
#define __FUNCT__ "readMatrixASCII"
PetscErrorCode readMatrixASCII(Mat& Matrix, char* fname, const PetscInt nrows, const PetscInt ncols, const char* type) {
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
PetscBool readInput(char* coordFname, char* phiFname, char* indicesFname, char* eqLengthFname, char* refFname,
		    char* oFname, PetscInt& numCoords, PetscInt& numCons, PetscInt& numCG)
{

    PetscFunctionBegin;
    PetscBool flg;

    PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "--ref", coordFname, PETSC_MAX_PATH_LEN, &flg);
    PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "--cg", phiFname, PETSC_MAX_PATH_LEN, &flg);
    PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "--indices", indicesFname, PETSC_MAX_PATH_LEN, &flg);
    PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "--lengths", eqLengthFname, PETSC_MAX_PATH_LEN, &flg);
    PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "--cgOP", refFname, PETSC_MAX_PATH_LEN, &flg);
    PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "--out", oFname, PETSC_MAX_PATH_LEN, &flg);

    numCoords = countLines(coordFname);
    numCG = countLines(phiFname);
    numCons = countLines(indicesFname);

    assert(numCons == countLines(eqLengthFname));

    PetscFunctionReturn(flg);
}

#undef __FUNCT__
#define __FUNCT__ "cleanUp"
PetscErrorCode cleanUp(Mat& Ref, Vec& EqLength, std::vector<Vec>& Coords, std::vector<Vec>& Phi ) {

    PetscFunctionBegin;
    PetscErrorCode ierr;

    ierr = MatDestroy(&Ref); CHKERRQ(ierr);
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
