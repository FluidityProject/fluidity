static char help[] = "2D Lineal Spring Analogy Smoother in serial and parallel.\n\n";
#include <petscksp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


void lin_smoother(int dimension, int num_nodes, int num_elements, int num_surf_elements, int num_owned_nodes, int * mapping, int * connectivity, double * phys_mesh, double * smooth_mesh, double * comp_mesh, int * surf_connectivity) {
  
  Vec            F,U_h;
  Mat            K;            
  KSP            ksp;
  PC             pc;
  PetscErrorCode ierr;
  PetscInt       i,n,m,its,Istart,Iend,Ii,len_new=0,j,local_surf_connectivity[dimension*num_surf_elements], global_surf_connectivity[dimension*num_surf_elements],
                 num_nodes_col_x[num_nodes],num_nodes_col_y[num_nodes];
  PetscScalar    length,x_disp[num_surf_elements],y_disp[num_surf_elements],z_disp[num_surf_elements],smoothed_x[num_nodes],smoothed_y[num_nodes],smoothed_z[num_nodes];

  double lij(int dimension, int nodei, int nodej){
    if (nodei == nodej){
    PetscPrintf(PETSC_COMM_WORLD,"Failure: nodei is the same as nodej.\n");
    length = 0.0;
    }
    else if (dimension == 1) {
      double xi = comp_mesh[dimension*nodei-1];
      double xj = comp_mesh[dimension*nodej-1];
      length = labs(xi-xj);
    }
    else if (dimension == 2) {
      double xi = comp_mesh[dimension*nodei]; double yi = comp_mesh[dimension*nodei+1];
      double xj = comp_mesh[dimension*nodej]; double yj = comp_mesh[dimension*nodej+1];
      length = sqrt(pow(xi-xj,2)+pow(yi-yj,2));
    }
    else {
      double xi = comp_mesh[dimension*nodei-3]; double yi = comp_mesh[dimension*nodei-2]; double zi = comp_mesh[dimension*nodei-1];
      double xj = comp_mesh[dimension*nodej-3]; double yj = comp_mesh[dimension*nodej-2]; double zj = comp_mesh[dimension*nodej-1];
      length = sqrt(pow(xi-xj,2)+pow(yi-yj,2)+pow(zi-zj,2));
    }
    
    return length;
  }

  /*Lineal Stiffness*/
  double * K_lin(int nodei, int nodej){
    nodej = nodej-1;
    static double K_lin_temp[16];
    double xi = comp_mesh[dimension*nodei]; double yi = comp_mesh[dimension*nodei+1];
    double xj = comp_mesh[dimension*nodej]; double yj = comp_mesh[dimension*nodej+1];
    double x_edge = xj - xi; double y_edge = yj - yi;
    double alpha = atan2(y_edge,x_edge);
    double inv_length = 1.0/lij(2,nodei,nodej);
    
    K_lin_temp[0]=inv_length*cos(alpha)*cos(alpha);K_lin_temp[1]=inv_length*sin(alpha)*cos(alpha);K_lin_temp[2]=-1.0*inv_length*cos(alpha)*cos(alpha);K_lin_temp[3]=-1.0*inv_length*cos(alpha)*sin(alpha);
    K_lin_temp[4]=inv_length*sin(alpha)*cos(alpha);K_lin_temp[5]=inv_length*sin(alpha)*sin(alpha);K_lin_temp[6]=-1.0*inv_length*sin(alpha)*cos(alpha);K_lin_temp[7]=-1.0*inv_length*sin(alpha)*sin(alpha);
    K_lin_temp[8]=-1.0*inv_length*cos(alpha)*cos(alpha);K_lin_temp[9]=-1.0*inv_length*sin(alpha)*cos(alpha);K_lin_temp[10]=inv_length*cos(alpha)*cos(alpha);K_lin_temp[11]=inv_length*sin(alpha)*cos(alpha);
    K_lin_temp[12]=-1.0*inv_length*sin(alpha)*cos(alpha);K_lin_temp[13]=-1.0*inv_length*sin(alpha)*sin(alpha);K_lin_temp[14]=inv_length*sin(alpha)*cos(alpha);K_lin_temp[15]=inv_length*sin(alpha)*sin(alpha);

    return K_lin_temp;
  }

  int conn_mat[num_elements][3];
  for(int m=0;m<num_elements;++m){
    conn_mat[m][0]=*(connectivity+(3*m));
    conn_mat[m][1]=*(connectivity+(3*m)+1);
    conn_mat[m][2]=*(connectivity+(3*m)+2);
  }


  /*Neighbour Matrix*/
  int neb_mat[num_owned_nodes][8];
  for(int nodei=0;nodei<num_owned_nodes;nodei++){
    int neb_a_count=0; int neb_b_count=0; int neb_c_count=0;
    int neb_a[8]={};int neb_b[8]={};int neb_c[8]={};int neb_tot[24]={};int neb_fin[8]={};
    
    for(int m=0;m<num_elements;++m){
     
      if (nodei+1 == conn_mat[m][0]){
	neb_a[neb_a_count]=conn_mat[m][1];
	neb_a[neb_a_count+1]=conn_mat[m][2];
	neb_a_count+=2;}
      
      if (nodei+1 == conn_mat[m][1]){
	neb_b[neb_b_count]=conn_mat[m][0];
	neb_b[neb_b_count+1]=conn_mat[m][2];
	neb_b_count+=2;}
      
      if (nodei+1 == conn_mat[m][2]){
	  neb_c[neb_c_count]=conn_mat[m][0];
	  neb_c[neb_c_count+1]=conn_mat[m][1];
	  neb_c_count+=2;}\
    }
    
    for(int n =0;n<8;n++){
      neb_tot[n]=neb_a[n];
      neb_tot[8+n]=neb_b[n];
      neb_tot[16+n]=neb_c[n];
    }

   //Removing duplicates 
   for(int m=0;m<24;++m){
     for(int n=m+1;n<24;++n){
       if(neb_tot[m]==neb_tot[n]){
    	  neb_tot[n]=0;
       }
     }
   }

   //sorting
   for(int i=1; i<=24-1; i++){
     for(int j=1; j<=24-i; j++){
       if(neb_tot[j-1] <= neb_tot[j]){
	 int t = neb_tot[j-1];
	 neb_tot[j-1]=neb_tot[j];
	 neb_tot[j]=t;
       }
     }
   }

   //Only interested in first 8 ints
   for(int m=0;m<8;m++){
     neb_mat[nodei][m]=neb_tot[m];}
  }

  
  MatCreate(PETSC_COMM_WORLD,&K);
  PetscObjectSetName((PetscObject) K, "Stiffness Matrix");
  MatSetSizes(K,2*num_owned_nodes,2*num_owned_nodes,PETSC_DETERMINE,PETSC_DETERMINE);
  MatSetUp(K);

  for(int Ii=0;Ii<num_owned_nodes;Ii++){
       for(int n=0;n<8;n++){
     	  int neb_hold = neb_mat[Ii][n];
	  if(neb_hold != 0){	    
	    MatSetValue(K,2*mapping[Ii],2*mapping[Ii],*(K_lin(Ii,neb_hold)+0),ADD_VALUES);
	    MatSetValue(K,2*mapping[Ii],2*mapping[Ii]+1,*(K_lin(Ii,neb_hold)+1),ADD_VALUES);
	    MatSetValue(K,2*mapping[Ii],2*mapping[neb_hold-1],*(K_lin(Ii,neb_hold)+2),ADD_VALUES);
	    MatSetValue(K,2*mapping[Ii],2*mapping[neb_hold-1]+1,*(K_lin(Ii,neb_hold)+3),ADD_VALUES);
	    
	    MatSetValue(K,2*mapping[Ii]+1,2*mapping[Ii],*(K_lin(Ii,neb_hold)+4),ADD_VALUES);
	    MatSetValue(K,2*mapping[Ii]+1,2*mapping[Ii]+1,*(K_lin(Ii,neb_hold)+5),ADD_VALUES);
	    MatSetValue(K,2*mapping[Ii]+1,2*mapping[neb_hold-1],*(K_lin(Ii,neb_hold)+6),ADD_VALUES);
	    MatSetValue(K,2*mapping[Ii]+1,2*mapping[neb_hold-1]+1,*(K_lin(Ii,neb_hold)+7),ADD_VALUES);

	    
	    MatSetValue(K,2*mapping[neb_hold-1],2*mapping[Ii],*(K_lin(Ii,neb_hold)+8),ADD_VALUES);
	    MatSetValue(K,2*mapping[neb_hold-1],2*mapping[Ii]+1,*(K_lin(Ii,neb_hold)+9),ADD_VALUES);
	    MatSetValue(K,2*mapping[neb_hold-1],2*mapping[neb_hold-1],*(K_lin(Ii,neb_hold)+10),ADD_VALUES);
	    MatSetValue(K,2*mapping[neb_hold-1],2*mapping[neb_hold-1]+1,*(K_lin(Ii,neb_hold)+11),ADD_VALUES);
	    
	    MatSetValue(K,2*mapping[neb_hold-1]+1,2*mapping[Ii],*(K_lin(Ii,neb_hold)+12),ADD_VALUES);
	    MatSetValue(K,2*mapping[neb_hold-1]+1,2*mapping[Ii]+1,*(K_lin(Ii,neb_hold)+13),ADD_VALUES);
	    MatSetValue(K,2*mapping[neb_hold-1]+1,2*mapping[neb_hold-1],*(K_lin(Ii,neb_hold)+14),ADD_VALUES);
	    MatSetValue(K,2*mapping[neb_hold-1]+1,2*mapping[neb_hold-1]+1,*(K_lin(Ii,neb_hold)+15),ADD_VALUES);
	  }
       }
  }
    
  MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);

  VecCreate(PETSC_COMM_WORLD,&F);
  PetscObjectSetName((PetscObject) F, "Load Vector");
  VecSetSizes(F,2*num_owned_nodes,PETSC_DECIDE);
  VecSetFromOptions(F);

  //Removing duplicate points in surf_connectivity
  for(int i=0;i<dimension*num_surf_elements;i++){
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

  
  for(int i=0;i<len_new;i++){
    global_surf_connectivity[2*i]=2*mapping[local_surf_connectivity[i]];
    global_surf_connectivity[2*i+1]=2*mapping[local_surf_connectivity[i]]+1;
  }
  

  MatZeroRows(K,2*len_new,global_surf_connectivity,1.0,NULL,NULL);

 for(int n=0;n<len_new;n++){
   x_disp[n]=phys_mesh[dimension*(local_surf_connectivity[n])] - comp_mesh[dimension*(local_surf_connectivity[n])];
   y_disp[n]=phys_mesh[dimension*(local_surf_connectivity[n])+1] - comp_mesh[dimension*(local_surf_connectivity[n])+1];
  }
 
 for(Ii=0;Ii<len_new;Ii++){
   if (local_surf_connectivity[Ii]<num_owned_nodes){
     VecSetValue(F,global_surf_connectivity[2*Ii],x_disp[Ii],INSERT_VALUES);
     VecSetValue(F,global_surf_connectivity[2*Ii+1],y_disp[Ii],INSERT_VALUES);
   }
 }
   
  VecAssemblyBegin(F);
  VecAssemblyEnd(F);

  MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);

  KSPCreate(PETSC_COMM_WORLD,&ksp);

  KSPSetOperators(ksp,K,K);

  KSPSetType(ksp,KSPGMRES);

  KSPGetPC(ksp,&pc);
  PCSetType(pc,PCSOR);
  KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
  KSPSetTolerances(ksp,1e-7,1e-7,PETSC_DEFAULT,PETSC_DEFAULT);
  KSPSetFromOptions(ksp);

  VecCreate(PETSC_COMM_WORLD,&U_h);
  PetscObjectSetName((PetscObject) U_h, "U_h");
  VecSetSizes(U_h,2*num_owned_nodes,PETSC_DECIDE);
  VecSetFromOptions(U_h);

  VecAssemblyBegin(U_h);
  VecAssemblyEnd(U_h);

  KSPSolve(ksp,F,U_h);
  KSPGetIterationNumber(ksp,&its);
  PetscPrintf(PETSC_COMM_WORLD,"ksp iter: %D\n", its);
  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp,&reason);
  PetscPrintf(PETSC_COMM_WORLD,"KSPConvergedReason: %D\n", reason);

  for(n=0;n<num_owned_nodes;n++){
     num_nodes_col_x[n]=2*mapping[n];
     num_nodes_col_y[n]=2*mapping[n]+1;
   }
  
  VecGetValues(U_h,num_owned_nodes,num_nodes_col_x,smoothed_x);
  VecGetValues(U_h,num_owned_nodes,num_nodes_col_y,smoothed_y);

  for(n=0;n<num_owned_nodes; n++){
    smooth_mesh[dimension*n+0] = smoothed_x[n]+comp_mesh[dimension*n+0];
    smooth_mesh[dimension*n+1] = smoothed_y[n]+comp_mesh[dimension*n+1];
   }

  
  VecDestroy(&F);VecDestroy(&U_h);MatDestroy(&K);KSPDestroy(&ksp);
}
