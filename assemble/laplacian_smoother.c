static char help[] = "Solves 1D, 2D and 3D Laplace's equation in serial and parallel.\n\n";
#include <petscksp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "solver_options.h"

void lap_smoother(int dimension, int num_nodes, int num_elements, int num_surf_elements, int num_owned_nodes, int* mapping, int * connectivity, double * phys_mesh, double* smooth_mesh, double * comp_mesh, int * surf_connectivity, Mat* K, struct solver_options* options, int debug_level) {
  
  Vec            F[3],U_h[3];          
  PetscInt       i,n,m,convert[4],*holder_new,
		 Ii,len_new=0,j,local_surf_connectivity[dimension*num_surf_elements], global_surf_connectivity[dimension*num_surf_elements];
  PetscScalar    p_value[4],A_holder[16],grad_holder[12],det,Fe,x_bound_points[num_surf_elements],Vol,A_inv_holder[4][4],length,
                 y_bound_points[num_surf_elements],z_bound_points[num_surf_elements],Area,Ke[4][4];

 //Decomposing F into F_x, F_y and F_z owing to nature of Laplace's equation (independence of variables) 
  VecCreate(PETSC_COMM_WORLD,&F[0]);
  PetscObjectSetName((PetscObject) F[0], "Fx");
  VecSetSizes(F[0],num_owned_nodes,PETSC_DECIDE);
  VecSetFromOptions(F[0]);
  if (dimension>1) {
    VecCreate(PETSC_COMM_WORLD,&F[1]);
    PetscObjectSetName((PetscObject) F[1], "Fy");
    VecSetSizes(F[1],num_owned_nodes,PETSC_DECIDE);
    VecSetFromOptions(F[1]);
    if (dimension>2) {
      VecCreate(PETSC_COMM_WORLD,&F[2]);
      PetscObjectSetName((PetscObject) F[2], "Fz");
      VecSetSizes(F[2],num_owned_nodes,PETSC_DECIDE);
      VecSetFromOptions(F[2]);
    }
  }


  //Assembly of Stiffness Matrix K
  if (dimension == 1){
    for(Ii=0;Ii<num_elements;Ii++){
      
      holder_new = &(connectivity[(dimension+1)*Ii]);
      convert[0]=holder_new[0]-1;
      convert[1]=holder_new[1]-1;
      int counter = 0;
      
    for(m=0;m<dimension+1;m++){
      p_value[0]=1.0;p_value[1]=comp_mesh[dimension*convert[m]];
      A_holder[counter]=p_value[0];A_holder[counter+1]=p_value[1];
      counter+=dimension+1;
     }
           
     det=A_holder[0]*A_holder[3] - A_holder[1]*A_holder[2];
 
     A_inv_holder[0][0]=(1/det)*A_holder[3];
     A_inv_holder[0][1]=-1.0*(1/det)*A_holder[1];
     A_inv_holder[1][0]=-1.0*(1/det)*A_holder[2];
     A_inv_holder[1][1]=(1/det)*A_holder[0];

     length = fabs(det);
  
     grad_holder[0]=A_inv_holder[1][0];grad_holder[1]=A_inv_holder[1][1];
          
     Ke[0][0]=length*(pow(grad_holder[0],2));
     Ke[0][1]=length*(grad_holder[0]*grad_holder[1]);
     Ke[1][0]=length*(grad_holder[1]*grad_holder[0]);
     Ke[1][1]=length*(pow(grad_holder[1],2));
     
     Fe=(length/2.0)*0.0;//As called for in Laplacian Smoothing

     if (holder_new[0]<=num_owned_nodes) {
       MatSetValue(*K,mapping[convert[0]],mapping[convert[0]],Ke[0][0],ADD_VALUES);
       MatSetValue(*K,mapping[convert[0]],mapping[convert[1]],Ke[0][1],ADD_VALUES);
       VecSetValue(F[0],mapping[convert[0]],Fe,ADD_VALUES);
     }
     if (holder_new[1]<=num_owned_nodes) {
       MatSetValue(*K,mapping[convert[1]],mapping[convert[0]],Ke[1][0],ADD_VALUES);
       MatSetValue(*K,mapping[convert[1]],mapping[convert[1]],Ke[1][1],ADD_VALUES);
       VecSetValue(F[0],mapping[convert[1]],Fe,ADD_VALUES);
     }
      
    }//End of 1D Assembly

  
  } else if (dimension == 2){//Beginning of 2D assembly
     for (Ii=0;Ii<num_elements;Ii++){

      holder_new = &(connectivity[(dimension+1)*Ii]);
      convert[0]=holder_new[0]-1;
      convert[1]=holder_new[1]-1;
      convert[2]=holder_new[2]-1;
      int counter = 0;
    
    for(m=0;m<dimension+1;m++){
      p_value[0]=1.0;p_value[1]=comp_mesh[dimension*convert[m]];p_value[2]=comp_mesh[dimension*convert[m]+1];
      A_holder[counter]=p_value[0];A_holder[counter+1]=p_value[1];A_holder[counter+2]=p_value[2];
      counter+=dimension+1;
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
       MatSetValue(*K,mapping[convert[0]],mapping[convert[0]],Ke[0][0],ADD_VALUES);
       MatSetValue(*K,mapping[convert[0]],mapping[convert[1]],Ke[0][1],ADD_VALUES);
       MatSetValue(*K,mapping[convert[0]],mapping[convert[2]],Ke[0][2],ADD_VALUES);
       VecSetValue(F[0],mapping[convert[0]],Fe,ADD_VALUES);
     }
     if (holder_new[1]<=num_owned_nodes) {
       MatSetValue(*K,mapping[convert[1]],mapping[convert[0]],Ke[1][0],ADD_VALUES);
       MatSetValue(*K,mapping[convert[1]],mapping[convert[1]],Ke[1][1],ADD_VALUES);
       MatSetValue(*K,mapping[convert[1]],mapping[convert[2]],Ke[1][2],ADD_VALUES);
       VecSetValue(F[0],mapping[convert[1]],Fe,ADD_VALUES);
     }
     if (holder_new[2]<=num_owned_nodes){
       MatSetValue(*K,mapping[convert[2]],mapping[convert[0]],Ke[2][0],ADD_VALUES);
       MatSetValue(*K,mapping[convert[2]],mapping[convert[1]],Ke[2][1],ADD_VALUES);
       MatSetValue(*K,mapping[convert[2]],mapping[convert[2]],Ke[2][2],ADD_VALUES);
       VecSetValue(F[0],mapping[convert[2]],Fe,ADD_VALUES);
     }
   } //End of 2D Assembly

  } else{//Beginning of 3D  
     for (Ii=0;Ii<num_elements;Ii++){

      holder_new = &(connectivity[(dimension+1)*Ii]);
      convert[0]=holder_new[0]-1;
      convert[1]=holder_new[1]-1;
      convert[2]=holder_new[2]-1;
      convert[3]=holder_new[3]-1;
      int counter = 0;

    for(m=0;m<dimension+1;m++){
      p_value[0]=1.0;p_value[1]=comp_mesh[dimension*convert[m]];p_value[2]=comp_mesh[dimension*convert[m]+1];p_value[3]=comp_mesh[dimension*convert[m]+2];
      A_holder[counter]=p_value[0];A_holder[counter+1]=p_value[1];A_holder[counter+2]=p_value[2];A_holder[counter+3]=p_value[3];
      counter+=dimension+1;
     }
//inverse  of a 4x4 matrix [Works thus far: shown to generate correct total volume of a unit cube]

    A_inv_holder[0][0] =  A_holder[5]*A_holder[10]*A_holder[15] - A_holder[5]*A_holder[11]*A_holder[14] - A_holder[9]*A_holder[6]*A_holder[15] + A_holder[9]*A_holder[7]*A_holder[14] + A_holder[13]*A_holder[6]*A_holder[11] - A_holder[13]*A_holder[7]*A_holder[10];

    A_inv_holder[0][1] = -A_holder[1]*A_holder[10]*A_holder[15] + A_holder[1]*A_holder[11]*A_holder[14] + A_holder[9]*A_holder[2]*A_holder[15] - A_holder[9]*A_holder[3]*A_holder[14] - A_holder[13]*A_holder[2]*A_holder[11] + A_holder[13]*A_holder[3]*A_holder[10];

    A_inv_holder[0][2] =  A_holder[1]*A_holder[6]*A_holder[15] - A_holder[1]*A_holder[7]*A_holder[14] - A_holder[5]*A_holder[2]*A_holder[15] + A_holder[5]*A_holder[3]*A_holder[14] + A_holder[13]*A_holder[2]*A_holder[7] - A_holder[13]*A_holder[3]*A_holder[6];

    A_inv_holder[0][3] = -A_holder[1]*A_holder[6]*A_holder[11] + A_holder[1]*A_holder[7]*A_holder[10] + A_holder[5]*A_holder[2]*A_holder[11] - A_holder[5]*A_holder[3]*A_holder[10] - A_holder[9]*A_holder[2]*A_holder[7] + A_holder[9]*A_holder[3]*A_holder[6];

    A_inv_holder[1][0] = -A_holder[4]*A_holder[10]*A_holder[15] + A_holder[4]*A_holder[11]*A_holder[14] + A_holder[8]*A_holder[6]*A_holder[15] - A_holder[8]*A_holder[7]*A_holder[14] - A_holder[12]*A_holder[6]*A_holder[11] + A_holder[12]*A_holder[7]*A_holder[10];

    A_inv_holder[1][1] =  A_holder[0]*A_holder[10]*A_holder[15] - A_holder[0]*A_holder[11]*A_holder[14] - A_holder[8]*A_holder[2]*A_holder[15] + A_holder[8]*A_holder[3]*A_holder[14] + A_holder[12]*A_holder[2]*A_holder[11] - A_holder[12]*A_holder[3]*A_holder[10];

    A_inv_holder[1][2] = -A_holder[0]*A_holder[6]*A_holder[15] + A_holder[0]*A_holder[7]*A_holder[14] + A_holder[4]*A_holder[2]*A_holder[15] - A_holder[4]*A_holder[3]*A_holder[14] - A_holder[12]*A_holder[2]*A_holder[7] + A_holder[12]*A_holder[3]*A_holder[6];

    A_inv_holder[1][3] =  A_holder[0]*A_holder[6]*A_holder[11] - A_holder[0]*A_holder[7]*A_holder[10] - A_holder[4]*A_holder[2]*A_holder[11] + A_holder[4]*A_holder[3]*A_holder[10] + A_holder[8]*A_holder[2]*A_holder[7] - A_holder[8]*A_holder[3]*A_holder[6];

    A_inv_holder[2][0] =  A_holder[4]*A_holder[9]*A_holder[15] - A_holder[4]*A_holder[11]*A_holder[13] - A_holder[8]*A_holder[5]*A_holder[15] + A_holder[8]*A_holder[7]*A_holder[13] + A_holder[12]*A_holder[5]*A_holder[11] - A_holder[12]*A_holder[7]*A_holder[9];

    A_inv_holder[2][1] = -A_holder[0]*A_holder[9]*A_holder[15] + A_holder[0]*A_holder[11]*A_holder[13] + A_holder[8]*A_holder[1]*A_holder[15] - A_holder[8]*A_holder[3]*A_holder[13] - A_holder[12]*A_holder[1]*A_holder[11] + A_holder[12]*A_holder[3]*A_holder[9];

    A_inv_holder[2][2] =  A_holder[0]*A_holder[5]*A_holder[15] - A_holder[0]*A_holder[7]*A_holder[13] - A_holder[4]*A_holder[1]*A_holder[15] + A_holder[4]*A_holder[3]*A_holder[13] + A_holder[12]*A_holder[1]*A_holder[7] - A_holder[12]*A_holder[3]*A_holder[5];

    A_inv_holder[2][3] = -A_holder[0]*A_holder[5]*A_holder[11] + A_holder[0]*A_holder[7]*A_holder[9] + A_holder[4]*A_holder[1]*A_holder[11] - A_holder[4]*A_holder[3]*A_holder[9] - A_holder[8]*A_holder[1]*A_holder[7] + A_holder[8]*A_holder[3]*A_holder[5];

    A_inv_holder[3][0] = -A_holder[4]*A_holder[9]*A_holder[14] + A_holder[4]*A_holder[10]*A_holder[13] + A_holder[8]*A_holder[5]*A_holder[14] - A_holder[8]*A_holder[6]*A_holder[13] - A_holder[12]*A_holder[5]*A_holder[10] + A_holder[12]*A_holder[6]*A_holder[9];

    A_inv_holder[3][1] =  A_holder[0]*A_holder[9]*A_holder[14] - A_holder[0]*A_holder[10]*A_holder[13] - A_holder[8]*A_holder[1]*A_holder[14] + A_holder[8]*A_holder[2]*A_holder[13] + A_holder[12]*A_holder[1]*A_holder[10] - A_holder[12]*A_holder[2]*A_holder[9];

    A_inv_holder[3][2] = -A_holder[0]*A_holder[5]*A_holder[14] + A_holder[0]*A_holder[6]*A_holder[13] + A_holder[4]*A_holder[1]*A_holder[14] - A_holder[4]*A_holder[2]*A_holder[13] - A_holder[12]*A_holder[1]*A_holder[6] + A_holder[12]*A_holder[2]*A_holder[5];

    A_inv_holder[3][3] =  A_holder[0]*A_holder[5]*A_holder[10] - A_holder[0]*A_holder[6]*A_holder[9] - A_holder[4]*A_holder[1]*A_holder[10] + A_holder[4]*A_holder[2]*A_holder[9] + A_holder[8]*A_holder[1]*A_holder[6] - A_holder[8]*A_holder[2]*A_holder[5];

    
   det = A_holder[0]*A_inv_holder[0][0] + A_holder[1]*A_inv_holder[1][0] + A_holder[2]*A_inv_holder[2][0] + A_holder[3]*A_inv_holder[3][0];
    
//NOTE: A_inv_holder = A_inv_holder * (1.f/det)
     for(i=0;i<4;i++){
          for(j=0;j<4;j++){
               A_inv_holder[i][j] = A_inv_holder[i][j] *(1.f/det);
                          }
   	             }	

     Vol = fabs(det)/6.0;
     
     grad_holder[0]=A_inv_holder[1][0];grad_holder[1]=A_inv_holder[1][1];grad_holder[2]=A_inv_holder[1][2];grad_holder[3]=A_inv_holder[1][3];
     grad_holder[4]=A_inv_holder[2][0];grad_holder[5]=A_inv_holder[2][1];grad_holder[6]=A_inv_holder[2][2];grad_holder[7]=A_inv_holder[2][3];
     grad_holder[8]=A_inv_holder[3][0];grad_holder[9]=A_inv_holder[3][1];grad_holder[10]=A_inv_holder[3][2];grad_holder[11]=A_inv_holder[3][3];

     //Ke = Vol*grad_holder'*grad_holder;
     
     Ke[0][0]=Vol*(pow(grad_holder[0],2)+pow(grad_holder[4],2)+pow(grad_holder[8],2));
     Ke[0][1]=Vol*(grad_holder[0]*grad_holder[1]+grad_holder[4]*grad_holder[5]+grad_holder[8]*grad_holder[9]);
     Ke[0][2]=Vol*(grad_holder[0]*grad_holder[2]+grad_holder[4]*grad_holder[6]+grad_holder[8]*grad_holder[10]);
     Ke[0][3]=Vol*(grad_holder[0]*grad_holder[3]+grad_holder[4]*grad_holder[7]+grad_holder[8]*grad_holder[11]);

     Ke[1][0]=Vol*(grad_holder[1]*grad_holder[0]+grad_holder[5]*grad_holder[4]+grad_holder[9]*grad_holder[8]);
     Ke[1][1]=Vol*(pow(grad_holder[1],2)+pow(grad_holder[5],2)+pow(grad_holder[9],2));
     Ke[1][2]=Vol*(grad_holder[1]*grad_holder[2]+grad_holder[5]*grad_holder[6]+grad_holder[9]*grad_holder[10]);
     Ke[1][3]=Vol*(grad_holder[1]*grad_holder[3]+grad_holder[5]*grad_holder[7]+grad_holder[9]*grad_holder[11]);

     Ke[2][0]=Vol*(grad_holder[2]*grad_holder[0]+grad_holder[6]*grad_holder[4]+grad_holder[10]*grad_holder[8]);
     Ke[2][1]=Vol*(grad_holder[2]*grad_holder[1]+grad_holder[6]*grad_holder[5]+grad_holder[10]*grad_holder[9]);
     Ke[2][2]=Vol*(pow(grad_holder[2],2)+pow(grad_holder[6],2)+pow(grad_holder[10],2));
     Ke[2][3]=Vol*(grad_holder[2]*grad_holder[3]+grad_holder[6]*grad_holder[7]+grad_holder[10]*grad_holder[11]);

     Ke[3][0]=Vol*(grad_holder[3]*grad_holder[0]+grad_holder[7]*grad_holder[4]+grad_holder[11]*grad_holder[8]);
     Ke[3][1]=Vol*(grad_holder[3]*grad_holder[1]+grad_holder[7]*grad_holder[5]+grad_holder[11]*grad_holder[9]);
     Ke[3][2]=Vol*(grad_holder[3]*grad_holder[2]+grad_holder[7]*grad_holder[6]+grad_holder[11]*grad_holder[10]);
     Ke[3][3]=Vol*(pow(grad_holder[3],2)+pow(grad_holder[7],2)+pow(grad_holder[11],2));

     Fe=0.0;//As called for in Laplacian Smoothing

     if (holder_new[0]<=num_owned_nodes){
       MatSetValue(*K,mapping[convert[0]],mapping[convert[0]],Ke[0][0],ADD_VALUES);
       MatSetValue(*K,mapping[convert[0]],mapping[convert[1]],Ke[0][1],ADD_VALUES);
       MatSetValue(*K,mapping[convert[0]],mapping[convert[2]],Ke[0][2],ADD_VALUES);
       MatSetValue(*K,mapping[convert[0]],mapping[convert[3]],Ke[0][3],ADD_VALUES);
       VecSetValue(F[0],mapping[convert[0]],Fe,ADD_VALUES);
     }
     if (holder_new[1]<=num_owned_nodes){
       MatSetValue(*K,mapping[convert[1]],mapping[convert[0]],Ke[1][0],ADD_VALUES);
       MatSetValue(*K,mapping[convert[1]],mapping[convert[1]],Ke[1][1],ADD_VALUES);
       MatSetValue(*K,mapping[convert[1]],mapping[convert[2]],Ke[1][2],ADD_VALUES);
       MatSetValue(*K,mapping[convert[1]],mapping[convert[3]],Ke[1][3],ADD_VALUES);
       VecSetValue(F[0],mapping[convert[1]],Fe,ADD_VALUES);
     }
     if (holder_new[2]<=num_owned_nodes){
       MatSetValue(*K,mapping[convert[2]],mapping[convert[0]],Ke[2][0],ADD_VALUES);
       MatSetValue(*K,mapping[convert[2]],mapping[convert[1]],Ke[2][1],ADD_VALUES);
       MatSetValue(*K,mapping[convert[2]],mapping[convert[2]],Ke[2][2],ADD_VALUES);
       MatSetValue(*K,mapping[convert[2]],mapping[convert[3]],Ke[2][3],ADD_VALUES);
       VecSetValue(F[0],mapping[convert[2]],Fe,ADD_VALUES);
     }
     if (holder_new[3]<=num_owned_nodes){
       MatSetValue(*K,mapping[convert[3]],mapping[convert[0]],Ke[3][0],ADD_VALUES);
       MatSetValue(*K,mapping[convert[3]],mapping[convert[1]],Ke[3][1],ADD_VALUES);
       MatSetValue(*K,mapping[convert[3]],mapping[convert[2]],Ke[3][2],ADD_VALUES);
       MatSetValue(*K,mapping[convert[3]],mapping[convert[3]],Ke[3][3],ADD_VALUES);
       VecSetValue(F[0],mapping[convert[3]],Fe,ADD_VALUES);
     }
   } //End of 3D assembly
  }
  MatAssemblyBegin(*K,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*K,MAT_FINAL_ASSEMBLY);

  //Removing duplicate points in surf_connectivity
  for(i=0; i< dimension*num_surf_elements; i++){

    if(surf_connectivity[i]>num_owned_nodes) continue;

    for(j=0;j< len_new;j++){

      if(surf_connectivity[i] == local_surf_connectivity[j]+1)
      break;
   }

     if (j==len_new)
        local_surf_connectivity[len_new++] = surf_connectivity[i]-1;
   }

  for(i=0;i<len_new;i++){
    global_surf_connectivity[i]=mapping[local_surf_connectivity[i]];
  }

  MatZeroRows(*K,len_new,global_surf_connectivity,1.0,NULL,NULL);
  
  for(n=0;n<len_new;n++){
    x_bound_points[n]=phys_mesh[dimension*(local_surf_connectivity[n])];
    if (dimension>1) {
      y_bound_points[n]=phys_mesh[dimension*(local_surf_connectivity[n])+1];
      if (dimension>2) {
	z_bound_points[n]=phys_mesh[dimension*(local_surf_connectivity[n])+2];
      }
    }
  }
  
  for(Ii=0;Ii<len_new;Ii++){
    if (local_surf_connectivity[Ii]<num_owned_nodes){
      VecSetValue(F[0],global_surf_connectivity[Ii],x_bound_points[Ii],INSERT_VALUES);
      if (dimension>1) {
	VecSetValue(F[1],global_surf_connectivity[Ii],y_bound_points[Ii],INSERT_VALUES);
	if (dimension>2) {
	  VecSetValue(F[2],global_surf_connectivity[Ii],z_bound_points[Ii],INSERT_VALUES);
	}
      }
    }
  }
    
  MatAssemblyBegin(*K,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*K,MAT_FINAL_ASSEMBLY);
  

  for (j=0; j<dimension; ++j){

    VecAssemblyBegin(F[j]);
    VecAssemblyEnd(F[j]);

    VecCreate(PETSC_COMM_WORLD,&(U_h[j]));
    PetscObjectSetName((PetscObject) U_h[j], "laplacian_smoother");
    VecSetSizes(U_h[j],num_owned_nodes,PETSC_DECIDE);
    VecSetFromOptions(U_h[j]);
    
    PetscScalar *aU_h;
    // Copy in first guess,
    VecGetArray(U_h[j],&aU_h);
    for(n=0;n<num_owned_nodes;n++){
      aU_h[n]=phys_mesh[dimension*n+j];
    }
    VecRestoreArray(U_h[j],&aU_h);
    
    VecAssemblyBegin(U_h[j]);
    VecAssemblyEnd(U_h[j]);
  }

  petsc_solve_many_c(U_h, *K, F, dimension, options, debug_level);

  for (j=0; j<dimension; ++j){
    // Copy back the result
    PetscScalar *aU_h;
    VecGetArray(U_h[j],&aU_h);
    for(n=0;n<num_owned_nodes;n++)
      {
	smooth_mesh[dimension*n+j] = aU_h[n];
      }
    VecRestoreArray(U_h[j],&aU_h);
    
    VecDestroy(&(F[j])); 
    VecDestroy(&(U_h[j]));    
  } 
}

