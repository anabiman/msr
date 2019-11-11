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

#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <fstream>
#include <string>
#include <petscksp.h>
#include <petscmat.h>
#include <vector>
#include <algorithm>

PetscErrorCode NewtonIter(std::vector<Vec>&, std::vector<Vec>&, PetscInt*, PetscInt*, Vec&, Mat&, Mat&, char*, PetscReal);
PetscErrorCode writeVector(std::vector<Vec>& Coords, const char* fname);

#define MAX_ITERS 20

#endif
