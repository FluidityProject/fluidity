static char help[] = "Solves 2D Laplacian equation in serial.\n\n";
#include <petscksp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void lap_smoother(int dimension, int num_nodes, int num_elements, int num_surf_elements, int * connectivity, double * phys_mesh, double* smooth_mesh, double * comp_mesh, int * surf_connectivity) {
  
  Vec            F,Fx,Fy,U_hx,U_hy,line_ele_vec;
  Mat            K,tri_ele_mat;            
  KSP            ksp_x,ksp_y; 
  PetscErrorCode ierr;
  PetscInt       i,n,m,col[3]={0,1,2},col_new[3],convert[3],line_total=num_surf_elements,c[num_surf_elements],
                 Istart,Iend,Ii,its,num_nodes_col[num_nodes],j_other=0,len_new=1,j;
  PetscMPIInt    size;
  PetscScalar    holder_new[3],p_value[3],A_holder[9],grad_holder[6],det,Fe,x_bound_points[num_surf_elements],
                 y_bound_points[num_surf_elements],Area,Ke[3][3],smoothed_x[num_nodes],smoothed_y[num_nodes],value[5];
  PetscViewer    viewer;
  float          A_inv_holder[3][3];

  
  //Creating PetscMat just for triangular elements
  MatCreate(PETSC_COMM_WORLD,&tri_ele_mat);
  PetscObjectSetName((PetscObject) tri_ele_mat, "triangle_connectivty");
  MatSetSizes(tri_ele_mat,PETSC_DECIDE,PETSC_DECIDE,num_elements,3);
  MatSetType(tri_ele_mat,MATAIJ);
  MatSetUp(tri_ele_mat);
  
  for (i=0;i<num_elements;i++){
    value[0]=connectivity[j_other]; value[1]=connectivity[j_other+1];value[2]=connectivity[j_other+2];
    MatSetValues(tri_ele_mat,1,&i,3,col,value,INSERT_VALUES);
    j_other+=3;
  }
  
  MatAssemblyBegin(tri_ele_mat,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(tri_ele_mat,MAT_FINAL_ASSEMBLY);
  
  MatCreate(PETSC_COMM_WORLD,&K);
  PetscObjectSetName((PetscObject) K, "Stiffness Matrix");
  MatSetSizes(K,PETSC_DECIDE,PETSC_DECIDE,num_nodes,num_nodes);
  MatSetUp(K);

  VecCreate(PETSC_COMM_WORLD,&F);
  PetscObjectSetName((PetscObject) F, "Load Vector");
  VecSetSizes(F,PETSC_DECIDE,num_nodes);
  VecSetFromOptions(F);

  MatGetOwnershipRange(tri_ele_mat,&Istart,&Iend);

  //Assembly of Stiffness Matrix K
  for (Ii=Istart;Ii<Iend;Ii++){

      MatGetValues(tri_ele_mat,1,&Ii,3,col,holder_new);
      convert[0]=holder_new[0]-1;convert[1]=holder_new[1]-1;convert[2]=holder_new[2]-1;
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

     MatSetValue(K,convert[0],convert[0],Ke[0][0],ADD_VALUES);MatSetValue(K,convert[0],convert[1],Ke[0][1],ADD_VALUES);
     MatSetValue(K,convert[0],convert[2],Ke[0][2],ADD_VALUES);MatSetValue(K,convert[1],convert[0],Ke[1][0],ADD_VALUES);
     MatSetValue(K,convert[1],convert[1],Ke[1][1],ADD_VALUES);MatSetValue(K,convert[1],convert[2],Ke[1][2],ADD_VALUES);
     MatSetValue(K,convert[2],convert[0],Ke[2][0],ADD_VALUES);MatSetValue(K,convert[2],convert[1],Ke[2][1],ADD_VALUES);
     MatSetValue(K,convert[2],convert[2],Ke[2][2],ADD_VALUES);
     
     VecSetValue(F,convert[0],Fe,ADD_VALUES);
     VecSetValue(F,convert[1],Fe,ADD_VALUES);
     VecSetValue(F,convert[2],Fe,ADD_VALUES);
    
   }
      
  MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(F);
  VecAssemblyEnd(F);

  //Removing duplicate points in in surf_connectivity
  for(i=1; i< dimension*num_surf_elements; i++){

   for(j=0; j< len_new ; j++)
   {

      if(surf_connectivity[i] == surf_connectivity[j])
      break;
   }

     if (j==len_new)
        surf_connectivity[len_new++] = surf_connectivity[i];
   }

 
  MatZeroRows(K,num_surf_elements,surf_connectivity-1,0.0,NULL,NULL);
  
  //Decomposing F into F_x and F_y, owing to nature of Laplace's equation (independence of variabls) 
  VecCreate(PETSC_COMM_WORLD,&Fx);
  PetscObjectSetName((PetscObject) Fx, "Fx");
  VecSetSizes(Fx,PETSC_DECIDE,num_nodes);
  VecSetFromOptions(Fx);
  VecDuplicate(Fx,&F);

  VecCreate(PETSC_COMM_WORLD,&Fy);
  PetscObjectSetName((PetscObject) Fy, "Fy");
  VecSetSizes(Fy,PETSC_DECIDE,num_nodes);
  VecSetFromOptions(Fy);
  VecDuplicate(Fy,&F);

  VecCreate(PETSC_COMM_WORLD,&line_ele_vec);
  PetscObjectSetName((PetscObject) line_ele_vec, "line_ele_vector");
  VecSetSizes(line_ele_vec,PETSC_DECIDE,num_surf_elements);
  VecSetFromOptions(line_ele_vec);
  VecGetOwnershipRange(line_ele_vec,&Istart,&Iend);
  
  for(n=0;n<num_surf_elements;n++){
    x_bound_points[n]=phys_mesh[2*(surf_connectivity[n]-1)];
    y_bound_points[n]=phys_mesh[(2*(surf_connectivity[n]-1))+1];
  }
  
  for(Ii=Istart;Ii<Iend;Ii++){

    MatSetValue(K,surf_connectivity[Ii]-1,surf_connectivity[Ii]-1,1.0,INSERT_VALUES);
    VecSetValue(Fx,surf_connectivity[Ii]-1,x_bound_points[Ii],INSERT_VALUES);
    VecSetValue(Fy,surf_connectivity[Ii]-1,y_bound_points[Ii],INSERT_VALUES);
  }

  
  MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);
  
  VecAssemblyBegin(Fx);
  VecAssemblyEnd(Fx);
  VecAssemblyBegin(Fy);
  VecAssemblyEnd(Fy);
  
  MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);

  KSPCreate(PETSC_COMM_WORLD,&ksp_x);
  KSPCreate(PETSC_COMM_WORLD,&ksp_y);

  KSPSetOperators(ksp_x,K,K);
  KSPSetOperators(ksp_y,K,K);

  KSPSetType(ksp_x,KSPGMRES);
  KSPSetType(ksp_y,KSPGMRES);
  
  KSPSetTolerances(ksp_x,1e-7,1e-7,PETSC_DEFAULT,PETSC_DEFAULT);
  KSPSetTolerances(ksp_y,1e-7,1e-7,PETSC_DEFAULT,PETSC_DEFAULT);
  
  KSPSetFromOptions(ksp_x);
  KSPSetFromOptions(ksp_y);

  VecCreate(PETSC_COMM_WORLD,&U_hx);
  PetscObjectSetName((PetscObject) U_hx, "U_hx");
  VecSetSizes(U_hx,PETSC_DECIDE,num_nodes);
  VecSetFromOptions(U_hx);
 
  VecCreate(PETSC_COMM_WORLD,&U_hy);
  PetscObjectSetName((PetscObject) U_hy, "U_hy");
  VecSetSizes(U_hy,PETSC_DECIDE,num_nodes);
  VecSetFromOptions(U_hy);
  
  KSPSolve(ksp_x,Fx,U_hx);
  KSPSolve(ksp_y,Fy,U_hy);
 
  KSPGetIterationNumber(ksp_x,&its);
  PetscPrintf(PETSC_COMM_WORLD,"ksp_x iter: %D\n", its);
  KSPGetIterationNumber(ksp_y,&its);
  PetscPrintf(PETSC_COMM_WORLD,"ksp_y iter: %D\n", its);

  
  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp_x,&reason);
  PetscPrintf(PETSC_COMM_WORLD,"KSPConvergedReasonX: %D\n", reason);
  KSPGetConvergedReason(ksp_y,&reason);
  PetscPrintf(PETSC_COMM_WORLD,"KSPConvergedReasonY: %D\n", reason);

  VecAssemblyBegin(U_hx);
  VecAssemblyEnd(U_hx);
  VecAssemblyBegin(U_hy);
  VecAssemblyEnd(U_hy);
  
  for(n=0;n<num_nodes;n++){
    num_nodes_col[n]=n;
  }
  VecGetValues(U_hx,num_nodes,num_nodes_col,smoothed_x);
  VecGetValues(U_hy,num_nodes,num_nodes_col,smoothed_y);

 for(n=0;n<num_nodes; n++){
    smooth_mesh[2*n+0] = smoothed_x[n];
    smooth_mesh[2*n+1] = smoothed_y[n];
  }

  VecDestroy(&Fx);VecDestroy(&Fy);VecDestroy(&U_hx);VecDestroy(&U_hy);VecDestroy(&line_ele_vec);
  MatDestroy(&K);MatDestroy(&tri_ele_mat);
  KSPDestroy(&ksp_x);KSPDestroy(&ksp_y);
 
}
