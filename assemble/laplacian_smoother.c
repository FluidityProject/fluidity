static char help[] = "Solves 2D Laplacian equation in serial.\n\n";
#include <petscksp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void lap_smoother(int dimension, int num_nodes, int num_elements, int num_surf_elements, int num_owned_nodes, int* mapping, int * connectivity, double * phys_mesh, double* smooth_mesh, double * comp_mesh, int * surf_connectivity) {
  
  Vec            F,Fx,Fy,U_hx,U_hy;
  Mat            K;            
  KSP            ksp_x;
  PC             pc;
  PetscErrorCode ierr;
  PetscInt       i,n,m,col[3]={0,1,2},col_new[3],convert[4],line_total=num_surf_elements,c[num_surf_elements],
		 Istart,Iend,Ii,its,num_nodes_col[num_nodes],j_other=0,len_new=0,j, local_surf_connectivity[dimension*num_surf_elements], global_surf_connectivity[dimension*num_surf_elements];
  PetscMPIInt    size;
  PetscInt       *holder_new;
  PetscScalar    p_value[4],A_holder[16],grad_holder[6],det,Fe,x_bound_points[num_surf_elements],Vol,
                 y_bound_points[num_surf_elements],Area,Ke[4][4],smoothed_x[num_nodes],smoothed_y[num_nodes],value[5];
  PetscViewer    viewer;
  float          A_inv_holder[4][4];
  
  MatCreate(PETSC_COMM_WORLD,&K);
  PetscObjectSetName((PetscObject) K, "Stiffness Matrix");
  MatSetSizes(K,num_owned_nodes,num_owned_nodes,PETSC_DETERMINE,PETSC_DETERMINE);
  MatSetUp(K);

  VecCreate(PETSC_COMM_WORLD,&F);
  PetscObjectSetName((PetscObject) F, "Load Vector");
  VecSetSizes(F,num_owned_nodes,PETSC_DECIDE);
  VecSetFromOptions(F);

  //Assembly of Stiffness Matrix K
  for (Ii=0;Ii<num_elements;Ii++){

    if (dimension == 1){

      holder_new = &(connectivity[(dimension+1)*Ii]);
      convert[0]=holder_new[0]-1;
      convert[1]=holder_new[1]-1;
      int counter = 0;

      //Do 1D FEM
    }

    else if (dimension == 2){

      holder_new = &(connectivity[(dimension+1)*Ii]);
      convert[0]=holder_new[0]-1;
      convert[1]=holder_new[1]-1;
      convert[2]=holder_new[2]-1;
      int counter = 0;
      
    for(m=0;m<3;m++){
      p_value[0]=1.0;p_value[1]=comp_mesh[2*convert[m]];p_value[2]=comp_mesh[2*convert[m]+1];
	A_holder[counter]=p_value[0];A_holder[counter+1]=p_value[1];A_holder[counter+2]=p_value[2];
	counter+=3;
     }
           
     det=A_holder[0]*(A_holder[4]*A_holder[8]-A_holder[5]*A_holder[7])-A_holder[1]*(A_holder[3]*A_holder[8]
		   - A_holder[5]*A_holder[6])+A_holder[2]*(A_holder[3]*A_holder[7]-A_holder[4]*A_holder[6]);
 
     A_inv_holder[0][0]=(1/det)*(A_holder[4]*A_holder[8]-A_holder[5]*A_holder[7]);
     A_inv_holder[0][1]=(1/det)*(A_holder[2]*A_holder[7]-A_holder[1]*A_holder[8]);A_inv_holder[0][2]=(1/det)*(A_holder[1]*A_holder[5]-A_holder[2]*A_holder[4]);
     A_inv_holder[1][0]=(1/det)*(A_holder[5]*A_holder[6]-A_holder[3]*A_holder[8]);
     A_inv_holder[1][1]=(1/det)*(A_holder[0]*A_holder[8]-A_holder[2]*A_holder[6]);A_inv_holder[1][2]=(1/det)*(A_holder[2]*A_holder[3]-A_holder[0]*A_holder[5]);
     A_inv_holder[2][0]=(1/det)*(A_holder[3]*A_holder[7]-A_holder[4]*A_holder[6]);
     A_inv_holder[2][1]=(1/det)*(A_holder[1]*A_holder[6]-A_holder[0]*A_holder[7]);A_inv_holder[2][2]=(1/det)*(A_holder[0]*A_holder[4]-A_holder[1]*A_holder[3]);

     Area = fabs(det)/2.0;
       
     grad_holder[0]=A_inv_holder[1][0];grad_holder[1]=A_inv_holder[1][1];grad_holder[2]=A_inv_holder[1][2];
     grad_holder[3]=A_inv_holder[2][0];grad_holder[4]=A_inv_holder[2][1];grad_holder[5]=A_inv_holder[2][2];
     
     Ke[0][0]=Area*(pow(grad_holder[0],2)+pow(grad_holder[3],2));Ke[0][1]=Area*(grad_holder[0]*grad_holder[1]+grad_holder[3]*grad_holder[4]);
     Ke[0][2]=Area*(grad_holder[0]*grad_holder[2]+grad_holder[3]*grad_holder[5]);Ke[1][0]=Area*(grad_holder[1]*grad_holder[0]+grad_holder[4]*grad_holder[3]);
     Ke[1][1]=Area*(pow(grad_holder[1],2)+pow(grad_holder[4],2));Ke[1][2]=Area*(grad_holder[1]*grad_holder[2]+grad_holder[4]*grad_holder[5]);
     Ke[2][0]=Area*(grad_holder[2]*grad_holder[0]+grad_holder[5]*grad_holder[3]);Ke[2][1]=Area*(grad_holder[2]*grad_holder[1]+grad_holder[5]*grad_holder[4]);
     Ke[2][2]=Area*(pow(grad_holder[2],2)+pow(grad_holder[5],2));

     Fe=(Area/3.0)*0.0;//As called for in Laplacian Smoothing


     if (holder_new[0]<=num_owned_nodes) {
       MatSetValue(K,mapping[convert[0]],mapping[convert[0]],Ke[0][0],ADD_VALUES);
       MatSetValue(K,mapping[convert[0]],mapping[convert[1]],Ke[0][1],ADD_VALUES);
       MatSetValue(K,mapping[convert[0]],mapping[convert[2]],Ke[0][2],ADD_VALUES);
       VecSetValue(F,mapping[convert[0]],Fe,ADD_VALUES);
     }
     if (holder_new[1]<=num_owned_nodes) {
       MatSetValue(K,mapping[convert[1]],mapping[convert[0]],Ke[1][0],ADD_VALUES);
       MatSetValue(K,mapping[convert[1]],mapping[convert[1]],Ke[1][1],ADD_VALUES);
       MatSetValue(K,mapping[convert[1]],mapping[convert[2]],Ke[1][2],ADD_VALUES);
       VecSetValue(F,mapping[convert[1]],Fe,ADD_VALUES);
     }
     if (holder_new[2]<=num_owned_nodes) {
       MatSetValue(K,mapping[convert[2]],mapping[convert[0]],Ke[2][0],ADD_VALUES);
       MatSetValue(K,mapping[convert[2]],mapping[convert[1]],Ke[2][1],ADD_VALUES);
       MatSetValue(K,mapping[convert[2]],mapping[convert[2]],Ke[2][2],ADD_VALUES);
       VecSetValue(F,mapping[convert[2]],Fe,ADD_VALUES);
     }
    }

    else {
     
      holder_new = &(connectivity[(dimension+1)*Ii]);
      convert[0]=holder_new[0]-1;
      convert[1]=holder_new[1]-1;
      convert[2]=holder_new[2]-1;
      convert[3]=holder_new[3]-1;
      int counter = 0;

    for(m=0;m<4;m++){
      p_value[0]=1.0;p_value[1]=comp_mesh[2*convert[m]];p_value[2]=comp_mesh[2*convert[m]+1];p_value[3]=comp_mesh[2*convert[m]+2];
      A_holder[counter]=p_value[0];A_holder[counter+1]=p_value[1];A_holder[counter+2]=p_value[2];A_holder[counter+3]=p_value[3];
      counter+=4;
     }
//determinant of a 4x4 matrix
      det = 
      A_holder[0]*(A_holder[5]*(A_holder[10]*A_holder[15]-A_holder[11]*A_holder[14]) -A_holder[6]*(A_holder[9]*A_holder[15]-A_holder[11]*A_holder[13]) + A_holder[7]*(A_holder[9]*A_holder[14]-A_holder[10]*A_holder[13]))
      - A_holder[1]*(A_holder[4]*(A_holder[10]*A_holder[15]-A_holder[11]*A_holder[14]) -A_holder[6]*(A_holder[8]*A_holder[15]-A_holder[11]*A_holder[12]) + A_holder[7]*(A_holder[8]*A_holder[14]-A_holder[10]*A_holder[12]))
      + A_holder[2]*(A_holder[4]*(A_holder[9]*A_holder[15]-A_holder[11]*A_holder[13]) -A_holder[5]*(A_holder[8]*A_holder[15]-A_holder[11]*A_holder[12]) + A_holder[7]*(A_holder[8]*A_holder[13]-A_holder[9]*A_holder[12]))
      - A_holder[3]*(A_holder[4]*(A_holder[9]*A_holder[14]-A_holder[10]*A_holder[13]) -A_holder[5]*(A_holder[8]*A_holder[14]-A_holder[10]*A_holder[12]) + A_holder[6]*(A_holder[8]*A_holder[13]- A_holder[9]*A_holder[12]));

      Vol = fabs(det)/6.0;

     if (holder_new[0]<=num_owned_nodes) {
       MatSetValue(K,mapping[convert[0]],mapping[convert[0]],Ke[0][0],ADD_VALUES);
       MatSetValue(K,mapping[convert[0]],mapping[convert[1]],Ke[0][1],ADD_VALUES);
       MatSetValue(K,mapping[convert[0]],mapping[convert[2]],Ke[0][2],ADD_VALUES);
       MatSetValue(K,mapping[convert[0]],mapping[convert[3]],Ke[0][3],ADD_VALUES);
       VecSetValue(F,mapping[convert[0]],Fe,ADD_VALUES);
     }
     if (holder_new[1]<=num_owned_nodes) {
       MatSetValue(K,mapping[convert[1]],mapping[convert[0]],Ke[1][0],ADD_VALUES);
       MatSetValue(K,mapping[convert[1]],mapping[convert[1]],Ke[1][1],ADD_VALUES);
       MatSetValue(K,mapping[convert[1]],mapping[convert[2]],Ke[1][2],ADD_VALUES);
       MatSetValue(K,mapping[convert[1]],mapping[convert[3]],Ke[1][3],ADD_VALUES);
       VecSetValue(F,mapping[convert[1]],Fe,ADD_VALUES);
     }
     if (holder_new[2]<=num_owned_nodes) {
       MatSetValue(K,mapping[convert[2]],mapping[convert[0]],Ke[2][0],ADD_VALUES);
       MatSetValue(K,mapping[convert[2]],mapping[convert[1]],Ke[2][1],ADD_VALUES);
       MatSetValue(K,mapping[convert[2]],mapping[convert[2]],Ke[2][2],ADD_VALUES);
       MatSetValue(K,mapping[convert[2]],mapping[convert[3]],Ke[2][3],ADD_VALUES);
       VecSetValue(F,mapping[convert[2]],Fe,ADD_VALUES);
     }
     if (holder_new[2]<=num_owned_nodes) {
       MatSetValue(K,mapping[convert[3]],mapping[convert[0]],Ke[3][0],ADD_VALUES);
       MatSetValue(K,mapping[convert[3]],mapping[convert[1]],Ke[3][1],ADD_VALUES);
       MatSetValue(K,mapping[convert[3]],mapping[convert[2]],Ke[3][2],ADD_VALUES);
       MatSetValue(K,mapping[convert[3]],mapping[convert[3]],Ke[3][3],ADD_VALUES);
       VecSetValue(F,mapping[convert[3]],Fe,ADD_VALUES);
     }
    }

   }
      
  MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(F);
  VecAssemblyEnd(F);
  //Removing duplicate points in in surf_connectivity
  for(i=0; i< dimension*num_surf_elements; i++){

    if(surf_connectivity[i]>num_owned_nodes) continue;

   for(j=0; j< len_new ; j++)
   {

      if(surf_connectivity[i] == local_surf_connectivity[j]+1)
      break;
   }

     if (j==len_new)
        local_surf_connectivity[len_new++] = surf_connectivity[i]-1;
   }

  for(i=0; i<len_new; i++) {
    global_surf_connectivity[i]=mapping[local_surf_connectivity[i]];
  }

  MatZeroRows(K,len_new,global_surf_connectivity,1.0,NULL,NULL);
  
  //Decomposing F into F_x and F_y, owing to nature of Laplace's equation (independence of variabls) 
  VecCreate(PETSC_COMM_WORLD,&Fx);
  PetscObjectSetName((PetscObject) Fx, "Fx");
  VecSetSizes(Fx,num_owned_nodes,PETSC_DECIDE);
  VecSetFromOptions(Fx);
  VecDuplicate(Fx,&F);

  VecCreate(PETSC_COMM_WORLD,&Fy);
  PetscObjectSetName((PetscObject) Fy, "Fy");
  VecSetSizes(Fy,num_owned_nodes,PETSC_DECIDE);
  VecSetFromOptions(Fy);
  VecDuplicate(Fy,&F);
  
  // VecCreate(PETSC_COMM_WORLD,&Fz);
  // PetscObjectSetName((PetscObject) Fz, "Fz");
  // VecSetSizes(Fz,num_owned_nodes,PETSC_DECIDE);
  // VecSetFromOptions(Fz);
  // VecDuplicate(Fz,&F);

  //3D boundary
  for(n=0;n<len_new;n++){
    x_bound_points[n]=phys_mesh[2*(local_surf_connectivity[n])];
    y_bound_points[n]=phys_mesh[2*(local_surf_connectivity[n])+1];
    // z_bound_points[n]=phys_mesh[2*(local_surf_connectivity[n])+2];
  }
  
  for(Ii=0;Ii<len_new;Ii++){
    if (local_surf_connectivity[Ii]<num_owned_nodes){
      VecSetValue(Fx,global_surf_connectivity[Ii],x_bound_points[Ii],INSERT_VALUES);
      VecSetValue(Fy,global_surf_connectivity[Ii],y_bound_points[Ii],INSERT_VALUES);
      // VecSetValue(Fz,global_surf_connectivity[Ii],z_bound_points[Ii],INSERT_VALUES);
   }
  }

  
  MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);
  
  VecAssemblyBegin(Fx);
  VecAssemblyEnd(Fx);
  VecAssemblyBegin(Fy);
  VecAssemblyEnd(Fy);
  // VecAssemblyBegin(Fz);
  // VecAssemblyEnd(Fz);

  KSPCreate(PETSC_COMM_WORLD,&ksp_x);

  KSPSetOperators(ksp_x,K,K);

  KSPSetType(ksp_x,KSPGMRES);

  KSPGetPC(ksp_x,&pc);
  PCSetType(pc,PCSOR);

  KSPSetInitialGuessNonzero(ksp_x,PETSC_TRUE);
  
  KSPSetTolerances(ksp_x,1e-7,1e-7,PETSC_DEFAULT,PETSC_DEFAULT);
  
  KSPSetFromOptions(ksp_x);

  VecCreate(PETSC_COMM_WORLD,&U_hx);
  PetscObjectSetName((PetscObject) U_hx, "U_hx");
  VecSetSizes(U_hx,num_owned_nodes,PETSC_DECIDE);
  VecSetFromOptions(U_hx);
 
  VecCreate(PETSC_COMM_WORLD,&U_hy);
  PetscObjectSetName((PetscObject) U_hy, "U_hy");
  VecSetSizes(U_hy,num_owned_nodes,PETSC_DECIDE);
  VecSetFromOptions(U_hy);

  // VecCreate(PETSC_COMM_WORLD,&U_hz);
  // PetscObjectSetName((PetscObject) U_hz, "U_hz");
  // VecSetSizes(U_hz,num_owned_nodes,PETSC_DECIDE);
  // VecSetFromOptions(U_hz);

  for(n=0;n<num_owned_nodes;n++){
    VecSetValue(U_hx,mapping[n],phys_mesh[2*n],INSERT_VALUES);
    VecSetValue(U_hy,mapping[n],phys_mesh[2*n+1],INSERT_VALUES);
    //  VecSetValue(U_hz,mapping[n],phys_mesh[2*n+2],INSERT_VALUES);

  }
  
  VecAssemblyBegin(U_hx);
  VecAssemblyEnd(U_hx);
  VecAssemblyBegin(U_hy);
  VecAssemblyEnd(U_hy);
  // VecAssemblyBegin(U_hz);
  // VecAssemblyEnd(U_hz);

  KSPSolve(ksp_x,Fx,U_hx);
  KSPGetIterationNumber(ksp_x,&its);
  PetscPrintf(PETSC_COMM_WORLD,"ksp_x iter: %D\n", its);
  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp_x,&reason);
  PetscPrintf(PETSC_COMM_WORLD,"KSPConvergedReasonX: %D\n", reason);

  KSPSolve(ksp_x,Fy,U_hy);
  KSPGetIterationNumber(ksp_x,&its);
  PetscPrintf(PETSC_COMM_WORLD,"ksp_y iter: %D\n", its);
  KSPGetConvergedReason(ksp_x,&reason);
  PetscPrintf(PETSC_COMM_WORLD,"KSPConvergedReasonY: %D\n", reason);

  // KSPSolve(ksp_x,Fz,U_hz);
  // KSPGetIterationNumber(ksp_x,&its);
  // PetscPrintf(PETSC_COMM_WORLD,"ksp_z iter: %D\n", its);
  // KSPGetConvergedReason(ksp_x,&reason);
  // PetscPrintf(PETSC_COMM_WORLD,"KSPConvergedReasonZ: %D\n", reason); 
  
  for(n=0;n<num_owned_nodes;n++){
    num_nodes_col[n]=mapping[n];
  }
  VecGetValues(U_hx,num_owned_nodes,num_nodes_col,smoothed_x);
  VecGetValues(U_hy,num_owned_nodes,num_nodes_col,smoothed_y);
  // VecGetValues(U_hz,num_owned_nodes,num_nodes_col,smoothed_z);

 for(n=0;n<num_owned_nodes; n++){
    smooth_mesh[2*n+0] = smoothed_x[n];
    smooth_mesh[2*n+1] = smoothed_y[n];
    // smooth_mesh[2*n+2] = smoothed_z[n];
  }

  VecDestroy(&Fx);VecDestroy(&Fy);//VecDestroy(&Fz);
  VecDestroy(&U_hx);VecDestroy(&U_hy);//VecDestroy(&U_hz);
  MatDestroy(&K);
  KSPDestroy(&ksp_x);
 
}

