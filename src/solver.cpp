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

/* This code does not run correctly in parallel. 
TODO: Find the bug in parallel mode.
*/

#undef __FUNCT__
#define __FUNCT__ "ComputeMetric"
PetscErrorCode ComputeMetric(std::vector<Vec> Coords, Mat& Adj, Vec& EqLength, Vec& Constraints) {
  /* ** Description **
     This function computes an array (distriuted on all procs) of equations satisfying
     the difference between the metric squared of two points in space and their
     equilibrium metric squared. This metric in this case corresponds to atomic
     distance (that is, for bonds/angles).

     |r_i - r_j|^2.0 - (l_ij)^2.0 with l_ij being the equilibrium metric to be satisfied.

     ** Credit **
     Andrew Abi Mansour
     Department of Chemistry
     Indiana University, Bloomington

     ** Last update **
     May 11, 2013 - 23:15
  */
	PetscFunctionBegin;

	PetscErrorCode ierr;

	// Extract appropriate dimension
	PetscInt numAtoms, Dims, numCons;
	Dims = Coords.size();
	VecGetSize(Coords[0], &numAtoms);
	VecGetSize(Constraints, &numCons);

	// Get local partition of the Constraints vector
	int istart, iend;
	VecGetOwnershipRange(Constraints,&istart,&iend);
	ierr = VecZeroEntries(Constraints);
	CHKERRQ(ierr);

	Vec consTmp;
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, numCons, &consTmp);

	for(auto dim = 0; dim < Dims; dim++) {
		
		ierr = MatMult(Adj, Coords[dim], consTmp); CHKERRQ(ierr);
		ierr = VecPow(consTmp, 2.0); CHKERRQ(ierr);

		VecAssemblyBegin(consTmp);
	        VecAssemblyEnd(consTmp);

		VecAXPY(Constraints, 1.0, consTmp);
		CHKERRQ(ierr);

	}

	ierr = VecAXPBY(Constraints, -1.0, 1.0, EqLength); // Cons = dx^y + dy^2 + dz^2 - Leq^2

	VecAssemblyBegin(Constraints);
	VecAssemblyEnd(Constraints);

	VecDestroy(&consTmp);

	PetscFunctionReturn(ierr);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeAtomicJacobi"
PetscErrorCode ComputeAtomicJacobi(std::vector<Vec>  Coords, Mat& Adj, PetscInt* indicesOne, PetscInt* indicesTwo, std::vector<Mat*> AtomicJacobi) {
  /* ** Description **
        This function constructs the Jacobian of the metric equations.

	** Size **
	The matrix should be of size NumCons x 3*NumAtoms

	** Last updated **

	** Credit **
	Andrew Abi Mansour
    Department of Chemistry
    Indiana University, Bloomington
  */
	PetscFunctionBegin;

	PetscErrorCode ierr;
	PetscInt istart, iend;
	PetscInt Dims, numAtoms, numCons;

	ierr = MatGetOwnershipRange(*AtomicJacobi[0], &istart, &iend); CHKERRQ(ierr);

	Dims = Coords.size();
	MatGetSize(Adj, &numCons, &numAtoms);

	PetscInt    *index            = new PetscInt[iend - istart];
  	PetscScalar *coordsDiffValues = new PetscScalar[iend - istart];

	Vec coordsDiff;
        ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, numCons, &coordsDiff);

	for(auto i = 0; i < iend - istart; i++)
		index[i] = i + istart;

	for(auto dim = 0; dim < Dims; dim++) {

                ierr = MatMult(Adj, Coords[dim], coordsDiff); CHKERRQ(ierr);
                ierr = VecScale(coordsDiff, 2.0); CHKERRQ(ierr);

                VecAssemblyBegin(coordsDiff);
                VecAssemblyEnd(coordsDiff);

		ierr = VecGetValues(coordsDiff, iend - istart, index, coordsDiffValues); CHKERRQ(ierr);

		for(auto i = istart; i < iend; i++) {

				ierr = MatSetValues(*AtomicJacobi[dim], 1, index + i - istart, 1, indicesOne + i - istart, coordsDiffValues + i - istart, INSERT_VALUES);
				CHKERRQ(ierr);

				PetscScalar valTmp = - coordsDiffValues[i- istart];
				ierr = MatSetValues(*AtomicJacobi[dim], 1, index + i - istart, 1, indicesTwo + i - istart, &valTmp, INSERT_VALUES);
				CHKERRQ(ierr);
			}

		ierr = MatAssemblyBegin(*AtomicJacobi[dim], MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  		ierr = MatAssemblyEnd(*AtomicJacobi[dim], MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	}

  	delete[] index;
	delete[] coordsDiffValues;
	VecDestroy(&coordsDiff);

	PetscFunctionReturn(ierr);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeAdjacencyMatrix"
PetscErrorCode ComputeAdjacencyMatrix(PetscInt* indicesOne, PetscInt* indicesTwo, Mat& Adj, const PetscInt Dims) {
  /*
     ** Description **
     This function computes an adjacency-like matrix which has the entries 1 and -1
     for each row: +1 for the first metric index, -1 for the second, i.e.

     \sum_d=1^3 Diag(P x r_d) x (P x r_d) = Metric^2.0

     ** Size **
     The total size of Adj should be NumCons x NumAtoms.

     ** Technical details **
     In the spirit of communication-free algorithms, Adj is assembled locally on each
     processor. For efficient storage and floating-point operations, the matrix is
     stored in sparse format.

     ** Last updated **
     May 10, 2013 - 11:45

     ** Credit **
     Andrew Abi Mansour
     Department of Chemistry
     Indiana University, Bloomington
  */
	PetscFunctionBegin;

	PetscErrorCode ierr;
	PetscInt NumCons, numAtoms, istart, iend;

	// Get the sizes of Indices (NumCons x 2)
	ierr = MatGetSize(Adj, &NumCons, &numAtoms);  CHKERRQ(ierr);
	ierr = MatGetOwnershipRange(Adj,&istart,&iend);  CHKERRQ(ierr);

	PetscInt    *index   = new PetscInt[iend - istart];
	PetscScalar *one     = new PetscScalar[iend - istart],
		    *none    = new PetscScalar[iend - istart];

	for(auto i = 0; i < iend - istart; i++) {

			index[i] = i + istart;
			one[i]   = 1.0;
			none[i]  = -1.0;
	}

	for(auto i = 0; i < iend - istart; i++) {
			ierr = MatSetValues(Adj, 1, index + i, 1, indicesOne + i, one + i, INSERT_VALUES); CHKERRQ(ierr);
			ierr = MatSetValues(Adj, 1, index + i, 1, indicesTwo + i, none + i, INSERT_VALUES); CHKERRQ(ierr);
		}

  	ierr = MatAssemblyBegin(Adj, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  	ierr = MatAssemblyEnd(Adj, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	//MatView(Adj, PETSC_VIEWER_STDOUT_WORLD);

  	delete[] index;
	delete[] one;
	delete[] none;

  	PetscFunctionReturn(ierr);
}

#undef __FUNCT__
#define __FUNCT__ "AssembleLagJacobi"
PetscErrorCode AssembleLagJacobi(std::vector<Vec> Coords, Mat& Adjacency, Mat& AdjacencyTrans, std::vector<Mat*> AtomicJacobi, Mat& LagJacobiGlobal) {
       PetscFunctionBegin;

       PetscErrorCode ierr;
       PetscInt Dims, numAtoms, NumCons, istart, iend;
       MatGetSize(AdjacencyTrans, &numAtoms, &NumCons);

       VecGetSize(Coords[0], &numAtoms);
       Dims = Coords.size();

       MatGetOwnershipRange(Adjacency, &istart, &iend);

       for(auto dim = 0; dim < Dims; dim++) {

	       // Compute Adj * Coords product
	       Vec Ar;
	       ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, NumCons, &Ar); CHKERRQ(ierr);
	       ierr = MatMult(Adjacency, Coords[dim], Ar); CHKERRQ(ierr);

	       ierr = VecScale(Ar, 2.0);

	       // Compute Jc * Adj.T product
	       Mat AdjAtom;
	       ierr = MatCreate(PETSC_COMM_WORLD, &AdjAtom); CHKERRQ(ierr);
	       ierr = MatSetSizes(AdjAtom, PETSC_DECIDE, PETSC_DECIDE, NumCons, NumCons); CHKERRQ(ierr);
	       ierr = MatSetType(AdjAtom, MATMPIAIJ); CHKERRQ(ierr);

	       ierr = MatMatMult(*AtomicJacobi[dim], AdjacencyTrans, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &AdjAtom); CHKERRQ(ierr);

	       ierr = MatDiagonalScale(AdjAtom, NULL, Ar); CHKERRQ(ierr);

	       if(dim == 0) {
		    		MatDuplicate(AdjAtom, MAT_COPY_VALUES, &LagJacobiGlobal);
	       		CHKERRQ(ierr);
	       }

	       else {
	            ierr = MatAXPY(LagJacobiGlobal, 1.0, AdjAtom, SAME_NONZERO_PATTERN);
	       	    CHKERRQ(ierr);
	       }

	       // Destroy allocated vectors & matrices
	       VecDestroy(&Ar);
	       MatDestroy(&AdjAtom);
	   }

       return ierr;
}

#undef __FUNCT__
#define __FUNCT__ "WriteMatrixToFile"
PetscErrorCode WriteMatrixToFile(Mat& Matrix, const char* ofname) {
      PetscFunctionBegin;
			PetscErrorCode ierr;

			PetscViewer viewer;
			ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF, ofname, &viewer); CHKERRQ(ierr);
			ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB); CHKERRQ(ierr);
			ierr = MatView(Matrix, viewer); CHKERRQ(ierr);
			ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

			PetscFunctionReturn(ierr);
}

#undef __FUNCT__
#define __FUNCT__ "NewtonIter"
PetscErrorCode NewtonIter(std::vector<Vec>& Coords, std::vector<Vec>& Phi, PetscInt* indicesOne, PetscInt* indicesTwo, Vec& EqLength, Mat& Ref, 
			  Mat& RefTransRef, char* oFname, PetscReal tol) {
     /*
        ** Description **
	This is the main computational engine that solves the equations

	\nabla_r f = r - r_u - Jc.T * mu + U.T * lambda = 0,
	\nabla_mu f = Diag(Adj * r) Adj * r - l_eq^2.0 = 0.
	\nabla_lambda f = \phi - U * r

	** Credit **
	Andrew Abi Mansour
	Department of Chemistry
	Indiana University, Bloomington
     */

	PetscFunctionBegin;

	PetscErrorCode ierr;
	PetscReal error1, error2;
	PetscInt NumCons, numAtoms, NumCG, Dims, ConsDims = 2, rank, size, istart, iend;

	MPI_Comm_size(PETSC_COMM_WORLD, &size);
  	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	Dims = Coords.size();
	VecGetSize(Coords[0], &numAtoms);
	VecGetSize(EqLength, &NumCons);
	VecGetSize(Phi[0], &NumCG);

  	Vec Constraints, mu, lambda, dr, CGCons;
 	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, NumCons, &Constraints); CHKERRQ(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, NumCG, &CGCons); CHKERRQ(ierr);
 	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, NumCons, &mu); CHKERRQ(ierr);
 	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, numAtoms, &dr); CHKERRQ(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, NumCG, &lambda); CHKERRQ(ierr);

	// Create the Seidel adjacency matrix locally (on each processor)
	Mat Adjacency;
	ierr = MatCreate(PETSC_COMM_WORLD, &Adjacency); CHKERRQ(ierr);
	ierr = MatSetSizes(Adjacency, PETSC_DECIDE, PETSC_DECIDE, NumCons, numAtoms); CHKERRQ(ierr);
	ierr = MatSetType(Adjacency, MATMPIAIJ); CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(Adjacency, ConsDims, PETSC_NULL, ConsDims, PETSC_NULL); CHKERRQ(ierr);
	ierr = ComputeAdjacencyMatrix(indicesOne, indicesTwo, Adjacency, Dims); CHKERRQ(ierr);

 	// Create AtomicJacobi and distribute it across all processors

	Mat AtomicJacobiX, AtomicJacobiY, AtomicJacobiZ;
	ierr = MatCreate(PETSC_COMM_WORLD, &AtomicJacobiX); CHKERRQ(ierr);
	ierr = MatSetSizes(AtomicJacobiX, PETSC_DECIDE, PETSC_DECIDE, NumCons, numAtoms); CHKERRQ(ierr);
	ierr = MatSetType(AtomicJacobiX, MATMPIAIJ); CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(AtomicJacobiX, ConsDims, PETSC_NULL, ConsDims, PETSC_NULL); CHKERRQ(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD, &AtomicJacobiY); CHKERRQ(ierr);
  	ierr = MatSetSizes(AtomicJacobiY, PETSC_DECIDE, PETSC_DECIDE, NumCons, numAtoms); CHKERRQ(ierr);
  	ierr = MatSetType(AtomicJacobiY, MATMPIAIJ); CHKERRQ(ierr);
  	ierr = MatMPIAIJSetPreallocation(AtomicJacobiY, ConsDims, PETSC_NULL, ConsDims, PETSC_NULL); CHKERRQ(ierr);

  	ierr = MatCreate(PETSC_COMM_WORLD, &AtomicJacobiZ); CHKERRQ(ierr);
  	ierr = MatSetSizes(AtomicJacobiZ, PETSC_DECIDE, PETSC_DECIDE, NumCons, numAtoms); CHKERRQ(ierr);
  	ierr = MatSetType(AtomicJacobiZ, MATMPIAIJ); CHKERRQ(ierr);
  	ierr = MatMPIAIJSetPreallocation(AtomicJacobiZ, ConsDims, PETSC_NULL, ConsDims, PETSC_NULL); CHKERRQ(ierr);

	std::vector<Mat*> AtomicJacobi(Dims);
	AtomicJacobi[0] = &AtomicJacobiX;
	AtomicJacobi[1]	= &AtomicJacobiY;
	AtomicJacobi[2]	= &AtomicJacobiZ;

	// perturb the coords vectors
	//Vec noise;
	//ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, numAtoms, &noise); CHKERRQ(ierr);

	//PetscRandom rctx;
	//ierr = PetscRandomCreate(PETSC_COMM_WORLD, &rctx); CHKERRQ(ierr);
	//ierr = PetscRandomSetType(rctx, PETSCRAND48);
	//ierr = VecSetRandom(noise, rctx); CHKERRQ(ierr);

	//ierr = VecScale(noise, 2.0); CHKERRQ(ierr);
	//ierr = VecShift(noise, -1.0); CHKERRQ(ierr);

	//for(auto it = Coords.begin(); it != Coords.end(); it++)
	//	VecAXPBY(*it, 1.0, 1.0, noise);

	//ierr = PetscRandomDestroy(&rctx); CHKERRQ(ierr);

	KSP ksp;
	//KSPSetTolerances(ksp,0.01,0.001,0.1,1000);
 	ierr = VecGetOwnershipRange(dr,&istart,&iend); CHKERRQ(ierr);

	Mat AdjacencyTrans;
	ierr = MatCreate(PETSC_COMM_WORLD, &AdjacencyTrans); CHKERRQ(ierr);
	ierr = MatSetSizes(AdjacencyTrans, PETSC_DECIDE, PETSC_DECIDE, numAtoms, NumCons); CHKERRQ(ierr);
	ierr = MatSetType(AdjacencyTrans, MATMPIAIJ); CHKERRQ(ierr);
	ierr = MatTranspose(Adjacency, MAT_INITIAL_MATRIX, &AdjacencyTrans); CHKERRQ(ierr);

	PetscScalar lagTmp;

	// Track prev coords of a newton iter
	std::vector<Vec> coordsRef(Dims);

	for(auto dim = 0; dim < Dims; dim++) {

			ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, numAtoms, &coordsRef[dim]);
			ierr = VecCopy(Coords[dim], coordsRef[dim]); CHKERRQ(ierr);

	}

	// SD parameters

	PetscScalar scaling = 1.0, sdScale = 0.1;
	PetscInt iter = 0, sdIters = 1;
	PetscScalar atomicDisp = tol;

	PetscScalar FSerror, CGerror, Lagrangian;
	
	try {
	 	while(atomicDisp >= tol) {

			//ierr = writeVector(Coords, oFname);

	 		ierr = ComputeMetric(Coords, Adjacency, EqLength, Constraints); CHKERRQ(ierr);

	 		ierr = ComputeAtomicJacobi(Coords, Adjacency, indicesOne, indicesTwo, AtomicJacobi); CHKERRQ(ierr);

	 		Mat LagJacobi;
			ierr = AssembleLagJacobi(Coords, Adjacency, AdjacencyTrans, AtomicJacobi, LagJacobi); CHKERRQ(ierr);
			ierr = MatShift(LagJacobi, scaling); CHKERRQ(ierr);

			// The nnz pattern does not change because the macromolecule topology does not change.

			// solve for mu
			// Create linear solver

			ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
			ierr = KSPSetInitialGuessNonzero(ksp, PETSC_FALSE); CHKERRQ(ierr);
			ierr = KSPSetType(ksp, KSPIBCGS); CHKERRQ(ierr);
			ierr = KSPSetFromOptions(ksp);

	 		ierr = KSPSetOperators(ksp, LagJacobi, LagJacobi); CHKERRQ(ierr);
	 		ierr = KSPSolve(ksp, Constraints, mu); CHKERRQ(ierr);
			ierr = KSPDestroy(&ksp); CHKERRQ(ierr);

			//VecView(Constraints, PETSC_VIEWER_STDOUT_SELF);

			// initialize values to monitor for analysis
			CGerror = .0;
			FSerror = .0;
			Lagrangian = .0;

	 		for(auto dim = 0; dim < Dims; dim++) {

					// Update r = ru + dr where dr = - Jt * mu
					ierr = MatMultTranspose(*AtomicJacobi[dim], mu, dr); CHKERRQ(ierr);
					ierr = VecAXPBY(Coords[dim], -1.0, 1.0, dr); CHKERRQ(ierr);
					VecNorm(dr, NORM_INFINITY, &atomicDisp);

					ierr = VecAssemblyBegin(Coords[dim]); CHKERRQ(ierr);
		  			ierr = VecAssemblyEnd(Coords[dim]); CHKERRQ(ierr);

					// solve for lambda vis SD
					
					for(auto sd = 0; sd < sdIters; sd++) {
					
						ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
						ierr = KSPSetInitialGuessNonzero(ksp, PETSC_FALSE); CHKERRQ(ierr);
						ierr = KSPSetType(ksp, KSPIBCGS); CHKERRQ(ierr);

						ierr = MatMultTranspose(Ref, Coords[dim], CGCons); CHKERRQ(ierr);

						ierr = VecAXPBY(CGCons, 1.0, -1.0, Phi[dim]); CHKERRQ(ierr);

						ierr = KSPSetOperators(ksp, RefTransRef, RefTransRef); CHKERRQ(ierr);
						ierr = KSPSetType(ksp, KSPIBCGS); CHKERRQ(ierr);
						ierr = KSPSolve(ksp, CGCons, lambda); CHKERRQ(ierr);
						ierr = KSPDestroy(&ksp); CHKERRQ(ierr);

						// Update r = ru + dr where dr = U * lambda
						ierr = MatMult(Ref, lambda, dr); CHKERRQ(ierr);
						ierr = VecAXPBY(Coords[dim], sdScale, 1.0, dr); CHKERRQ(ierr);

		 				ierr = VecAssemblyBegin(Coords[dim]); CHKERRQ(ierr);
		  				ierr = VecAssemblyEnd(Coords[dim]); CHKERRQ(ierr);

				 		VecNorm(dr, NORM_INFINITY, &error1); 
						atomicDisp += error1 * sdScale;

					}

					ierr = VecDot(lambda, CGCons, &lagTmp); CHKERRQ(ierr);
					Lagrangian += lagTmp;
					
					ierr = VecAXPY(coordsRef[dim], -1.0, Coords[dim]); CHKERRQ(ierr);
					ierr = VecDot(coordsRef[dim], coordsRef[dim], &lagTmp); CHKERRQ(ierr);
                                        Lagrangian += 0.5 * lagTmp;

					ierr = VecCopy(Coords[dim], coordsRef[dim]); CHKERRQ(ierr);

					VecNorm(Constraints, NORM_INFINITY, &error2);
					FSerror += error2;

					VecNorm(CGCons, NORM_INFINITY, &error2);
					CGerror += error2;
	  		}

				ierr = VecDot(mu, Constraints, &lagTmp); CHKERRQ(ierr);
				Lagrangian += lagTmp;

				if(!rank) {

					std::cout <<  "atom disp = " << atomicDisp << ", CG error = " << CGerror / 3.0 << " , FS error = " << FSerror / 3.0 << std::endl;
	
					std::ofstream myfile;
					myfile.open ("error.dat", std::ios::app);
  					myfile << CGerror / 3.0 << " ," << FSerror / 3.0 << " , " << Lagrangian << std::endl;
  					myfile.close();
				}

				MatDestroy(&LagJacobi);
				iter++;

				if(iter >MAX_ITERS)
					break;
	 		}
		}
		catch(const std::exception& e) {
			std::cout << "Failure in function " << __FUNCT__ << " with error " << e.what() << std::endl;
		}

 	//MatView(LagJacobi, PETSC_VIEWER_STDOUT_WORLD);
  	//VecView(Constraints,PETSC_VIEWER_STDOUT_SELF);

	VecNorm(CGCons, NORM_INFINITY, &error1);
        VecNorm(Constraints, NORM_INFINITY, &error2);

	if(!rank)
        	std::cout << "Scheme converged in " << iter << " iterations with CG error = " << error1 << " , and Micro error = " << error2 << std::endl;

	for(auto dim = 0; dim < Dims; dim++) {
			ierr = MatDestroy(AtomicJacobi[dim]); CHKERRQ(ierr);
	}

	ierr = MatDestroy(&AdjacencyTrans); CHKERRQ(ierr);
	ierr = MatDestroy(&Adjacency); CHKERRQ(ierr);

 	ierr = VecDestroy(&Constraints); CHKERRQ(ierr);
	ierr = VecDestroy(&dr); CHKERRQ(ierr);
	ierr = VecDestroy(&mu); CHKERRQ(ierr);

	CHKERRQ(ierr);

	PetscFunctionReturn(ierr);
}

