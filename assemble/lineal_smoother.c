static char help[] = "2D Lineal Spring Analogy Smoother in serial and parallel.\n\n";
#include <petscksp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void lin_smoother(int dimension, int num_nodes, int num_elements, int num_surf_elements, int num_owned_nodes, int * mapping, int * connectivity, double * phys_mesh, double * smooth_mesh, double * comp_mesh, int * surf_connectivity) {
  
  Vec            F,Fx,Fy,Fz,U_hx,U_hy,U_hz;
  Mat            K;            
  KSP            ksp_x;
  PC             pc;
  PetscErrorCode ierr;
  PetscInt       i,n,m,col[3]={0,1,2},col_new[3],convert[4],line_total=num_surf_elements,c[num_surf_elements],*holder_new,
		 Istart,Iend,Ii,its,num_nodes_col[num_nodes],j_other=0,len_new=0,j,local_surf_connectivity[dimension*num_surf_elements], global_surf_connectivity[dimension*num_surf_elements];
  PetscScalar    p_value[4],A_holder[16],grad_holder[12],det,Fe,x_bound_points[num_surf_elements],Vol,A_inv_holder[4][4],length,
                 y_bound_points[num_surf_elements],z_bound_points[num_surf_elements],Area,Ke[4][4],smoothed_x[num_nodes],smoothed_y[num_nodes],smoothed_z[num_nodes],value[5];

  MatCreate(PETSC_COMM_WORLD,&K);
  PetscObjectSetName((PetscObject) K, "Stiffness Matrix");
  MatSetSizes(K,2*num_owned_nodes,2*num_owned_nodes,PETSC_DETERMINE,PETSC_DETERMINE);
  MatSetUp(K);

  VecCreate(PETSC_COMM_WORLD,&F);
  PetscObjectSetName((PetscObject) F, "Load Vector");
  VecSetSizes(F,2*num_owned_nodes,PETSC_DECIDE);
  VecSetFromOptions(F);
  
  
  MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(F);
  VecAssemblyEnd(F);
 
   VecDestroy(&F);
   MatDestroy(&K);
 }
