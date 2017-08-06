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

#include <petscksp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "solver_options.h"

void setup_ksp_from_options_c(KSP ksp, struct solver_options *options) {

// Initialise a PETSc KSP solver from a set of solver options 
//
// This routine duplicates much of the functionality of the equivalent 
// Fortran subroutine 

  PC pc;

  KSPSetType(ksp,options->ksptype);
  KSPGetPC(ksp,&pc);

  int size;
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  if (strcmp(options->pctype,PCHYPRE)) { 
    PCSetType(pc,options->pctype);
    PCHYPRESetType(pc, options->hypretype);} 
  else if (strcmp(options->pctype,PCSOR) 
	   || strcmp(options->pctype,PCEISENSTAT)) {
    if (size == 1) {
      PCSetType(pc,options->pctype);
    } else {
      PCSetType(pc, PCBJACOBI);
      PCSetUp(pc);
      KSP *subksps;
      PC subpc;
      PCBJacobiGetSubKSP(pc, NULL, NULL, &subksps);
      KSPGetPC(subksps[0], &subpc);
      PCSetType(subpc,options->pctype);
    }
  }
  KSPSetInitialGuessNonzero(ksp,~options->start_from_zero);
  KSPSetTolerances(ksp,options->rtol,options->atol,PETSC_DEFAULT,options->max_its);
  KSPSetFromOptions(ksp);
}

void petsc_solve_many_c(Vec *x, Mat M, Vec *b, int neqns,
		   struct solver_options *options,
		   int debug_level) {
  KSP ksp;
 int n, its;
  KSPCreate(PETSC_COMM_WORLD,&ksp);

  KSPSetOperators(ksp,M,M);

  setup_ksp_from_options_c(ksp, options);
  
  char profiler_title[80];
  const char     *matname;
  PetscObjectGetName((PetscObject)(x[0]),&matname);
  strncpy(profiler_title,matname,80);
  strncat(profiler_title,"::solve",80);

  cprofiler_tic(profiler_title);
  for (n=0;n<neqns;++n) {
    KSPSolve(ksp, b[n], x[n]);
  }
  cprofiler_toc(profiler_title);
  KSPGetIterationNumber(ksp, &its);
  if (debug_level>1) PetscPrintf(PETSC_COMM_WORLD,"ksp iter: %D\n", its);
  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp, &reason);
  if (debug_level>1) PetscPrintf(PETSC_COMM_WORLD,"KSPConvergedReason: %D\n", reason);

  KSPDestroy(&ksp);
}

void petsc_solve_c(Vec x, Mat M, Vec b,
		   struct solver_options *options,
		   int debug_level) {
  petsc_solve_many_c(&x,M,&b,1,options, debug_level);
}
