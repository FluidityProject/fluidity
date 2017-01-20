static char help[] = "2D/3D Lineal Spring Analogy Smoother in serial and parallel.\n\n";
#include <petscksp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


void lin_smoother(int dimension, int num_nodes, int num_elements, int num_surf_elements, int num_owned_nodes, int * mapping, int connectivity[][3], double * phys_mesh, double * smooth_mesh, double * comp_mesh, int * surf_connectivity, int* findrm, int* colm, Mat* K) {
  
  Vec            F,U_h;
  KSP            ksp;
  PC             pc;
  PetscErrorCode ierr;
  PetscInt       i,n,m,j,its,nodei,Ii,len_new=0,local_surf_connectivity[dimension*num_surf_elements], global_surf_connectivity[dimension*num_surf_elements],
    num_nodes_col_x[num_nodes],num_nodes_col_y[num_nodes],num_nodes_col_z[num_nodes];
  PetscScalar    length,x_disp[num_surf_elements],y_disp[num_surf_elements],z_disp[num_surf_elements],smoothed_x[num_nodes],smoothed_y[num_nodes],smoothed_z[num_nodes];

  double lij(int nodei, int nodej){
    double xi,yi,zi,xj,yj,zj,length;
    if (nodei == nodej){
    PetscPrintf(PETSC_COMM_WORLD,"Failure: nodei is the same as nodej.\n");
    length = 0.0;
    }
    else if (dimension == 1) {
      xi = comp_mesh[dimension*nodei-1];
      xj = comp_mesh[dimension*nodej-1];
      length = labs(xi-xj);
    }
    else if (dimension == 2) {
      xi = comp_mesh[dimension*nodei]; yi = comp_mesh[dimension*nodei+1];
      xj = comp_mesh[dimension*nodej]; yj = comp_mesh[dimension*nodej+1];
      length = sqrt(pow(xi-xj,2)+pow(yi-yj,2));
    }
    else {
      xi = comp_mesh[dimension*nodei];yi = comp_mesh[dimension*nodei+1];zi = comp_mesh[dimension*nodei+2];
      xj = comp_mesh[dimension*nodej];yj = comp_mesh[dimension*nodej+1];zj = comp_mesh[dimension*nodej+2];
      length = sqrt(pow(xi-xj,2)+pow(yi-yj,2)+pow(zi-zj,2));
    }
    
    return length;
  }

  /*Lineal Stiffness*/
  double * K_lin(int nodei, int nodej){
    
    nodej = nodej-1;
    int n;
    double xi,yi,xj,yj,x_edge,y_edge,alpha,inv_length;

    xi = comp_mesh[dimension*nodei];yi = comp_mesh[dimension*nodei+1];
    xj = comp_mesh[dimension*nodej];yj = comp_mesh[dimension*nodej+1];
    x_edge = xj - xi;y_edge = yj - yi;
    alpha = atan2(y_edge,x_edge);inv_length = 1.0/lij(nodei,nodej);
    
    if(dimension == 2){
      static double K_lin_temp[3];
      K_lin_temp[0]=-1.0*inv_length*cos(alpha)*cos(alpha);
      K_lin_temp[1]=-1.0*inv_length*cos(alpha)*sin(alpha);
      K_lin_temp[2]=-1.0*inv_length*sin(alpha)*sin(alpha);

      return K_lin_temp;
    }
    else{
      double zi,zj,z_edge,phi;
      zi=comp_mesh[dimension*nodei+2];zj=comp_mesh[dimension*nodei+2];z_edge=zj-zi;phi=acos(z_edge/lij(nodei,nodej));
      static double K_lin_temp_3d[8];


      K_lin_temp_3d[0]=-1.0*inv_length*cos(alpha)*cos(alpha)*sin(phi);
      K_lin_temp_3d[1]=-1.0*inv_length*sin(alpha)*cos(alpha)*sin(phi);
      K_lin_temp_3d[2]=-1.0*inv_length*cos(alpha)*cos(phi); 
      K_lin_temp_3d[3]=-1.0*inv_length*sin(alpha)*sin(alpha)*sin(phi);
      K_lin_temp_3d[4]=-1.0*inv_length*sin(alpha)*cos(phi);
      K_lin_temp_3d[5]=-1.0*inv_length*cos(alpha)*sin(phi)*cos(phi);
      K_lin_temp_3d[6]=-1.0*inv_length*sin(alpha)*sin(phi)*cos(phi);
      K_lin_temp_3d[7]=-1.0*inv_length*cos(phi)*cos(phi);

      return K_lin_temp_3d;

    }
  }

  for(Ii=0;Ii<num_owned_nodes;Ii++){
    
    if(dimension==2){
  
       for(n=findrm[Ii];n<findrm[Ii+1];++n){
	 
	 int neb_hold = colm[n-1];
	 double * tmp = K_lin(Ii,neb_hold);
	 
	 MatSetValue(*K,2*mapping[Ii],2*mapping[neb_hold-1],*(tmp+0),ADD_VALUES);
	 MatSetValue(*K,2*mapping[Ii],2*mapping[neb_hold-1]+1,*(tmp+1),ADD_VALUES);
	 MatSetValue(*K,2*mapping[Ii]+1,2*mapping[neb_hold-1],*(tmp+1),ADD_VALUES);
	 MatSetValue(*K,2*mapping[Ii]+1,2*mapping[neb_hold-1]+1,*(tmp+2),ADD_VALUES);

	 MatSetValue(*K,2*mapping[neb_hold-1],2*mapping[neb_hold-1],-1.0*(*(tmp+0)),ADD_VALUES);
	 MatSetValue(*K,2*mapping[neb_hold-1],2*mapping[neb_hold-1]+1,-1.0*(*(tmp+1)),ADD_VALUES);
	 MatSetValue(*K,2*mapping[neb_hold-1]+1,2*mapping[neb_hold-1],-1.0*(*(tmp+1)),ADD_VALUES);
	 MatSetValue(*K,2*mapping[neb_hold-1]+1,2*mapping[neb_hold-1]+1,-1.0*(*(tmp+2)),ADD_VALUES);
	    
       }
    }
    else{
       for(n=findrm[Ii];n<findrm[Ii+1];++n){
	 
	 int neb_hold = colm[n-1];
	 double * tmp = K_lin(Ii,neb_hold);
	 
	 MatSetValue(*K,3*mapping[Ii],3*mapping[neb_hold-1],*(tmp+0),ADD_VALUES);
	 MatSetValue(*K,3*mapping[Ii],3*mapping[neb_hold-1]+1,*(tmp+1),ADD_VALUES);
	 MatSetValue(*K,3*mapping[Ii],3*mapping[neb_hold-1]+2,*(tmp+2),ADD_VALUES);
	 MatSetValue(*K,3*mapping[Ii]+1,3*mapping[neb_hold-1],*(tmp+1),ADD_VALUES);
	 MatSetValue(*K,3*mapping[Ii]+1,3*mapping[neb_hold-1]+1,*(tmp+3),ADD_VALUES);
	 MatSetValue(*K,3*mapping[Ii]+1,3*mapping[neb_hold-1]+2,*(tmp+4),ADD_VALUES);

	 MatSetValue(*K,3*mapping[Ii]+2,3*mapping[neb_hold-1],*(tmp+5),ADD_VALUES);
	 MatSetValue(*K,3*mapping[Ii]+2,3*mapping[neb_hold-1]+1,*(tmp+6),ADD_VALUES);
	 MatSetValue(*K,3*mapping[Ii]+2,3*mapping[neb_hold-1]+2,*(tmp+7),ADD_VALUES);
	 MatSetValue(*K,3*mapping[neb_hold-1],3*mapping[neb_hold-1],-1.0*(*(tmp+0)),ADD_VALUES);
	 MatSetValue(*K,3*mapping[neb_hold-1],3*mapping[neb_hold-1]+1,-1.0*(*(tmp+1)),ADD_VALUES);
	 MatSetValue(*K,3*mapping[neb_hold-1],3*mapping[neb_hold-1]+2,-1.0*(*(tmp+2)),ADD_VALUES);

	 MatSetValue(*K,3*mapping[neb_hold-1]+1,3*mapping[neb_hold-1],-1.0*(*(tmp+1)),ADD_VALUES);
	 MatSetValue(*K,3*mapping[neb_hold-1]+1,3*mapping[neb_hold-1]+1,-1.0*(*(tmp+3)),ADD_VALUES);
	 MatSetValue(*K,3*mapping[neb_hold-1]+1,3*mapping[neb_hold-1]+2,-1.0*(*(tmp+4)),ADD_VALUES);
	 MatSetValue(*K,3*mapping[neb_hold-1]+2,3*mapping[neb_hold-1],-1.0*(*(tmp+5)),ADD_VALUES);
	 MatSetValue(*K,3*mapping[neb_hold-1]+2,3*mapping[neb_hold-1]+1,-1.0*(*(tmp+6)),ADD_VALUES);
	 MatSetValue(*K,3*mapping[neb_hold-1]+2,3*mapping[neb_hold-1]+2,-1.0*(*(tmp+7)),ADD_VALUES);

       }
    }
  }
    
  MatAssemblyBegin(*K,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*K,MAT_FINAL_ASSEMBLY);

  VecCreate(PETSC_COMM_WORLD,&F);
  PetscObjectSetName((PetscObject) F, "Load Vector");
  VecSetSizes(F,dimension*num_owned_nodes,PETSC_DECIDE);
  VecSetFromOptions(F);

  //Removing duplicate points in surf_connectivity
  for(i=0;i<dimension*num_surf_elements;i++){
    if(surf_connectivity[i]>num_owned_nodes){
      continue;}
    
    for(j=0;j<len_new;j++)
      {
	if(surf_connectivity[i] == local_surf_connectivity[j]+1){
	  break;}
      }
    if(j==len_new){
      local_surf_connectivity[len_new++] = surf_connectivity[i]-1;
    }
  }

  if(dimension==2){
    for(i=0;i<len_new;i++){
      global_surf_connectivity[dimension*i]=dimension*mapping[local_surf_connectivity[i]];
      global_surf_connectivity[dimension*i+1]=dimension*mapping[local_surf_connectivity[i]]+1;
    }
  }
  else{
    for(i=0;i<len_new;i++){
      global_surf_connectivity[dimension*i]=dimension*mapping[local_surf_connectivity[i]];
      global_surf_connectivity[dimension*i+1]=dimension*mapping[local_surf_connectivity[i]]+1;
      global_surf_connectivity[dimension*i+2]=dimension*mapping[local_surf_connectivity[i]]+2;
    }
  }

  MatZeroRows(*K,dimension*len_new,global_surf_connectivity,1.0,NULL,NULL);

 if(dimension==2){
 for(n=0;n<len_new;n++){
   x_disp[n]=phys_mesh[dimension*(local_surf_connectivity[n])] - comp_mesh[dimension*(local_surf_connectivity[n])];
   y_disp[n]=phys_mesh[dimension*(local_surf_connectivity[n])+1] - comp_mesh[dimension*(local_surf_connectivity[n])+1];
  }
 
 for(Ii=0;Ii<len_new;Ii++){
   if (local_surf_connectivity[Ii]<num_owned_nodes){
     VecSetValue(F,global_surf_connectivity[dimension*Ii],x_disp[Ii],INSERT_VALUES);
     VecSetValue(F,global_surf_connectivity[dimension*Ii+1],y_disp[Ii],INSERT_VALUES);
   }
 }
  }
  else{
    for(n=0;n<len_new;n++){
      x_disp[n]=phys_mesh[dimension*(local_surf_connectivity[n])] - comp_mesh[dimension*(local_surf_connectivity[n])];
      y_disp[n]=phys_mesh[dimension*(local_surf_connectivity[n])+1] - comp_mesh[dimension*(local_surf_connectivity[n])+1];
      z_disp[n]=phys_mesh[dimension*(local_surf_connectivity[n])+2] - comp_mesh[dimension*(local_surf_connectivity[n])+2];
    }
    for(Ii=0;Ii<len_new;Ii++){
      if (local_surf_connectivity[Ii]<num_owned_nodes){
	VecSetValue(F,global_surf_connectivity[dimension*Ii],x_disp[Ii],INSERT_VALUES);
	VecSetValue(F,global_surf_connectivity[dimension*Ii+1],y_disp[Ii],INSERT_VALUES);
	VecSetValue(F,global_surf_connectivity[dimension*Ii+2],z_disp[Ii],INSERT_VALUES);
      }
    }
  }
  VecAssemblyBegin(F);
  VecAssemblyEnd(F);

  MatAssemblyBegin(*K,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*K,MAT_FINAL_ASSEMBLY);

  KSPCreate(PETSC_COMM_WORLD,&ksp);

  KSPSetOperators(ksp,*K,*K);

  KSPSetType(ksp,KSPGMRES);

  KSPGetPC(ksp,&pc);
  PCSetType(pc,PCSOR);
  KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
  KSPSetTolerances(ksp,1e-7,1e-7,PETSC_DEFAULT,PETSC_DEFAULT);
  KSPSetFromOptions(ksp);

  VecCreate(PETSC_COMM_WORLD,&U_h);
  PetscObjectSetName((PetscObject) U_h, "U_h");
  VecSetSizes(U_h,dimension*num_owned_nodes,PETSC_DECIDE);
  VecSetFromOptions(U_h);

  VecAssemblyBegin(U_h);
  VecAssemblyEnd(U_h);

  KSPSolve(ksp,F,U_h);
  KSPGetIterationNumber(ksp,&its);
  PetscPrintf(PETSC_COMM_WORLD,"ksp iter: %D\n", its);
  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp,&reason);
  PetscPrintf(PETSC_COMM_WORLD,"KSPConvergedReason: %D\n", reason);

  if(dimension==2){
    for(n=0;n<num_owned_nodes;n++){
      num_nodes_col_x[n]=dimension*mapping[n];
      num_nodes_col_y[n]=dimension*mapping[n]+1;
    }
    VecGetValues(U_h,num_owned_nodes,num_nodes_col_x,smoothed_x);
    VecGetValues(U_h,num_owned_nodes,num_nodes_col_y,smoothed_y);

    for(n=0;n<num_owned_nodes; n++){
      smooth_mesh[dimension*n+0] = smoothed_x[n]+comp_mesh[dimension*n+0];
      smooth_mesh[dimension*n+1] = smoothed_y[n]+comp_mesh[dimension*n+1];
    }
  }
  else{
    for(n=0;n<num_owned_nodes;n++){
      num_nodes_col_x[n]=dimension*mapping[n];
      num_nodes_col_y[n]=dimension*mapping[n]+1;
      num_nodes_col_z[n]=dimension*mapping[n]+2;
    }
    VecGetValues(U_h,num_owned_nodes,num_nodes_col_x,smoothed_x);
    VecGetValues(U_h,num_owned_nodes,num_nodes_col_y,smoothed_y);
    VecGetValues(U_h,num_owned_nodes,num_nodes_col_z,smoothed_z);

    for(n=0;n<num_owned_nodes;n++){
      smooth_mesh[dimension*n+0] = smoothed_x[n]+comp_mesh[dimension*n+0];
      smooth_mesh[dimension*n+1] = smoothed_y[n]+comp_mesh[dimension*n+1];
      smooth_mesh[dimension*n+2] = smoothed_z[n]+comp_mesh[dimension*n+2];
    }

  }

  
  VecDestroy(&F);VecDestroy(&U_h);KSPDestroy(&ksp);
}
