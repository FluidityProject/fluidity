/* !    Copyright (C) 2006 Imperial College London and others. */
/* !     */
/* !    Please see the AUTHORS file in the main source directory for a full list */
/* !    of copyright holders. */
/* ! */
/* !    Prof. C Pain */
/* !    Applied Modelling and Computation Group */
/* !    Department of Earth Science and Engineering */
/* !    Imperial College London */
/* ! */
/* !    amcgsoftware@imperial.ac.uk */
/* !     */
/* !    This library is free software; you can redistribute it and/or */
/* !    modify it under the terms of the GNU Lesser General Public */
/* !    License as published by the Free Software Foundation, */
/* !    version 2.1 of the License. */
/* ! */
/* !    This library is distributed in the hope that it will be useful, */
/* !    but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* !    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/* !    Lesser General Public License for more details. */
/* ! */
/* !    You should have received a copy of the GNU Lesser General Public */
/* !    License along with this library; if not, write to the Free Software */
/* !    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 */
/* !    USA */

#ifndef SOLVER_OPTIONS_H
#define SOLVER_OPTIONS_H
#include <petscksp.h>

struct solver_options {
  char ksptype[80];
  PetscInt restart;
  char pctype[80];
  char hypretype[80];
  PetscReal  atol;
  PetscReal rtol;
  PetscInt max_its;
  PetscBool start_from_zero;
};

void setup_ksp_from_options_c(KSP ksp, struct solver_options *options);
void petsc_solve_many_c(Vec *x, Mat M, Vec *b, int neqns,
		   struct solver_options *options,
		   int debug_level);
void petsc_solve_c(Vec x, Mat M, Vec b,
		   struct solver_options *options,
		   int debug_level);

#endif
