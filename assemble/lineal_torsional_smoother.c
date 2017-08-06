static char help[] = "2D Lineal Torsional Spring Analogy Smoother in serial and parallel.\n\n";
#include <petscksp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "solver_options.h"


void lin_tor_smoother(int dimension, int num_nodes, int num_elements, int num_surf_elements, int num_owned_nodes, int * mapping, int * connectivity, double * phys_mesh, double * smooth_mesh, double * comp_mesh, int * surf_connectivity, int* findrm, int* colm, Mat* K, struct solver_options* options, int debug_level) {
  
  Vec            F,U_h;         
  PetscInt       i,n,j,Ii,len_new=0,local_surf_connectivity[dimension*num_surf_elements], global_surf_connectivity[dimension*num_surf_elements];
  PetscScalar    x_disp[num_surf_elements],y_disp[num_surf_elements],z_disp[num_surf_elements];

   MatSetOption(*K,MAT_SYMMETRIC,PETSC_TRUE);
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

  double alpha_ang(int nodei, int nodej){
    double alpha;
    double xi = comp_mesh[dimension*nodei]; double yi = comp_mesh[dimension*nodei+1];
    double xj = comp_mesh[dimension*nodej]; double yj = comp_mesh[dimension*nodej+1];
    double x_edge = xj - xi; double y_edge = yj - yi;
    alpha = atan2(y_edge,x_edge);
    return alpha;
  }
    
  /*Lineal Stiffness*/
  double * K_lin(int nodei, int nodej){
    double xi,yi,xj,yj,x_edge,y_edge,alpha,inv_length;

    xi = comp_mesh[dimension*nodei];yi = comp_mesh[dimension*nodei+1];
    xj = comp_mesh[dimension*nodej];yj = comp_mesh[dimension*nodej+1];
    x_edge = xj - xi;y_edge = yj - yi;
    alpha = atan2(y_edge,x_edge);inv_length = 1.0/lij(nodei,nodej);
    
    static double K_lin_temp[3];
    K_lin_temp[0]=-1.0*inv_length*cos(alpha)*cos(alpha);
    K_lin_temp[1]=-1.0*inv_length*cos(alpha)*sin(alpha);
    K_lin_temp[2]=-1.0*inv_length*sin(alpha)*sin(alpha);
    
    return K_lin_temp;
  }

  double area(int nodei, int nodej, int nodek){
    double s = (lij(nodei,nodej) + lij(nodei,nodek) + lij(nodej,nodek))/2.0;
    double area = sqrt(s*(s-lij(nodei,nodej))*(s-lij(nodei,nodek))*(s-lij(nodej,nodek)));
    return area;
  }

  double tors(int nodei, int nodej, int nodek){
    double tors = (pow(lij(nodei,nodej),2) * pow(lij(nodei,nodek),2))/(4*pow(area(nodei,nodej,nodek),2));
    return tors;
  }
  
  double * C(int nodei, int nodej, int nodek){
    double* C_temp;
    C_temp = malloc(sizeof(double)*9);
    
    C_temp[0]=tors(nodei,nodej,nodek);
    C_temp[1]=0;C_temp[2]=0;C_temp[3]=0;
    C_temp[4]=tors(nodej,nodek,nodei);
    C_temp[5]=0;C_temp[6]=0;C_temp[7]=0;
    C_temp[8]=tors(nodek,nodei,nodej);
    return C_temp;
  }

  double x_edge(int nodei,int nodej){
    double xi,xj,x_edge;
    xi = comp_mesh[dimension*nodei];
    xj = comp_mesh[dimension*nodej];
    x_edge = xj - xi;
    return x_edge;
  }

  double y_edge(int nodei,int nodej){
    double yi,yj,y_edge;
    yi = comp_mesh[dimension*nodei+1];
    yj = comp_mesh[dimension*nodej+1];
    y_edge = yj - yi;
    return y_edge;
  }

  double a_rot(int nodei, int nodej){
    double a_rot = x_edge(nodei,nodej)/(pow(lij(nodei,nodej),2));
    return a_rot;
  }

  double b_rot(int nodei, int nodej){
    double b_rot = y_edge(nodei,nodej)/(pow(lij(nodei,nodej),2));
    return b_rot;
  }

  double * rot(int nodei, int nodej, int nodek){
    double* rot_mat;
    rot_mat = malloc(sizeof(double)*18);

    rot_mat[0]=b_rot(nodei,nodek)-b_rot(nodei,nodej);
    rot_mat[1]=a_rot(nodei,nodej)-a_rot(nodei,nodek);
    rot_mat[2]=b_rot(nodei,nodej);
    rot_mat[3]=-1.0*a_rot(nodei,nodej);
    rot_mat[4]=-1.0*b_rot(nodei,nodek);
    rot_mat[5]=a_rot(nodei,nodek);
    
    rot_mat[6]=-1.0*b_rot(nodej,nodei);
    rot_mat[7]=a_rot(nodej,nodei);
    rot_mat[8]=b_rot(nodej,nodei)-b_rot(nodej,nodek);
    rot_mat[9]=a_rot(nodej,nodek)-a_rot(nodej,nodei);
    rot_mat[10]=b_rot(nodej,nodek);
    rot_mat[11]=-1.0*a_rot(nodej,nodek);

    rot_mat[12]=b_rot(nodek,nodei);
    rot_mat[13]=-1.0*a_rot(nodek,nodei);
    rot_mat[14]=-1.0*b_rot(nodek,nodej);
    rot_mat[15]=a_rot(nodek,nodej);
    rot_mat[16]=b_rot(nodek,nodej)-b_rot(nodek,nodei);
    rot_mat[17]=a_rot(nodek,nodei)-a_rot(nodek,nodej);
    
    return rot_mat;
  }
  
  double * K_tor(int nodei, int nodej, int nodek){
    int i,j,k,r1=6,c1=3,c2=3;
    double rot_mat[3][6],rot_trans_mat[6][3], C_mat[3][3];
    double *A_mat,*K_tor_mat,*rot_hold, *C_hold,*K_tor_mat_trim;
    rot_hold = rot(nodei,nodej,nodek);
    C_hold = C(nodei,nodej,nodek);
    
    K_tor_mat = malloc(sizeof(double)*36);
    K_tor_mat_trim = malloc(sizeof(double)*24);
    A_mat = malloc(sizeof(double)*18);

    rot_mat[0][0]=rot_hold[0]; rot_mat[0][1]=rot_hold[1]; rot_mat[0][2]=rot_hold[2]; rot_mat[0][3]=rot_hold[3]; rot_mat[0][4]=rot_hold[4]; rot_mat[0][5]=rot_hold[5];
    rot_mat[1][0]=rot_hold[6]; rot_mat[1][1]=rot_hold[7]; rot_mat[1][2]=rot_hold[8]; rot_mat[1][3]=rot_hold[9]; rot_mat[1][4]=rot_hold[10];rot_mat[1][5]=rot_hold[11];
    rot_mat[2][0]=rot_hold[12];rot_mat[2][1]=rot_hold[13];rot_mat[2][2]=rot_hold[14];rot_mat[2][3]=rot_hold[15];rot_mat[2][4]=rot_hold[16];rot_mat[2][5]=rot_hold[17];
      
    rot_trans_mat[0][0]=rot_hold[0]; rot_trans_mat[0][1]=rot_hold[6]; rot_trans_mat[0][2]=rot_hold[12];
    rot_trans_mat[1][0]=rot_hold[1]; rot_trans_mat[1][1]=rot_hold[7]; rot_trans_mat[1][2]=rot_hold[13];
    rot_trans_mat[2][0]=rot_hold[2]; rot_trans_mat[2][1]=rot_hold[8]; rot_trans_mat[2][2]=rot_hold[14];
    rot_trans_mat[3][0]=rot_hold[3]; rot_trans_mat[3][1]=rot_hold[9]; rot_trans_mat[3][2]=rot_hold[15];
    rot_trans_mat[4][0]=rot_hold[4]; rot_trans_mat[4][1]=rot_hold[10];rot_trans_mat[4][2]=rot_hold[16];
    rot_trans_mat[5][0]=rot_hold[5]; rot_trans_mat[5][1]=rot_hold[11];rot_trans_mat[5][2]=rot_hold[17];

    C_mat[0][0]=C_hold[0]; C_mat[0][1]=C_hold[1]; C_mat[0][2]=C_hold[2];
    C_mat[1][0]=C_hold[3]; C_mat[1][1]=C_hold[4]; C_mat[1][2]=C_hold[5];
    C_mat[2][0]=C_hold[6]; C_mat[2][1]=C_hold[7]; C_mat[2][2]=C_hold[8];

    for(i=0;i<r1;i++){
      for(j=0;j<c2;j++){
	A_mat[(3*i+j)]=0;
      }
    }
    for(i=0;i<r1;i++){
      for(j=0;j<c2;j++){
	for(k=0;k<c1;k++){
	  A_mat[(3*i+j)]+=rot_trans_mat[i][k]*C_mat[k][j];
	}
      }
    }
    
    r1=6;c1=3;c2=6;
    for(i=0;i<r1;i++){
      for(j=0;j<c2;j++){
	K_tor_mat[(6*i+j)]=0;
      }
    }
    
    for(i=0;i<r1;i++){
      for(j=0;j<c2;j++){
	for(k=0;k<c1;k++){
	  K_tor_mat[(6*i+j)]+=A_mat[(3*i+k)]*rot_mat[k][j];
	}
      }
    }
    free(A_mat);free(C_hold);free(rot_hold);

    for(i=0;i<4;i++){
      for(j=0;j<6;j++){
	K_tor_mat_trim[(6*i+j)]=K_tor_mat[(6*i+j)];
      }
    }
    free(K_tor_mat);
    return K_tor_mat_trim;
  }
  
  for(Ii=0;Ii<num_owned_nodes;Ii++){
    for(n=findrm[Ii];n<findrm[Ii+1];++n){
      int neb_hold = colm[n-1]-1;
      double * tmp = K_lin(Ii,neb_hold);

      MatSetValue(*K,mapping[Ii],mapping[neb_hold],*(tmp+0),ADD_VALUES);
      MatSetValue(*K,mapping[Ii],mapping[num_nodes+neb_hold],*(tmp+1),ADD_VALUES);
      MatSetValue(*K,mapping[num_nodes+Ii],mapping[neb_hold],*(tmp+1),ADD_VALUES);
      MatSetValue(*K,mapping[num_nodes+Ii],mapping[num_nodes+neb_hold],*(tmp+2),ADD_VALUES);
      
      tmp = K_lin(neb_hold,Ii);
      
      MatSetValue(*K,mapping[Ii],mapping[Ii],-1.0*(*(tmp+0)),ADD_VALUES);
      MatSetValue(*K,mapping[Ii],mapping[num_nodes+Ii],-1.0*(*(tmp+1)),ADD_VALUES);
      MatSetValue(*K,mapping[num_nodes+Ii],mapping[Ii],-1.0*(*(tmp+1)),ADD_VALUES);
      MatSetValue(*K,mapping[num_nodes+Ii],mapping[num_nodes+Ii],-1.0*(*(tmp+2)),ADD_VALUES);
    }	    
  }

  int ele ;
  for(ele=0;ele<num_elements;++ele){
    int face;
    for (face=0;face<3;++face) {

      Ii=connectivity[3*ele+face]-1;
      if (Ii>=num_owned_nodes){continue;}
      int neb_hold = connectivity[3*ele+(face+1)%3]-1;
      int k1 = connectivity[3*ele+(face+2)%3]-1;

      if (Ii<num_owned_nodes) {

	double *K_tor1_holder=K_tor(Ii,neb_hold,k1);	    

      	MatSetValue(*K,mapping[Ii],mapping[Ii],(*(K_tor1_holder+0)),ADD_VALUES);
      	MatSetValue(*K,mapping[Ii],mapping[num_nodes+Ii],(*(K_tor1_holder+1)),ADD_VALUES);
      	MatSetValue(*K,mapping[Ii],mapping[neb_hold],(*(K_tor1_holder+2)),ADD_VALUES);
      	MatSetValue(*K,mapping[Ii],mapping[num_nodes+neb_hold],(*(K_tor1_holder+3)),ADD_VALUES);
      	MatSetValue(*K,mapping[Ii],mapping[k1],(*(K_tor1_holder+4)),ADD_VALUES);
      	MatSetValue(*K,mapping[Ii],mapping[num_nodes+k1],(*(K_tor1_holder+5)),ADD_VALUES);

      	MatSetValue(*K,mapping[num_nodes+Ii],mapping[Ii],(*(K_tor1_holder+6)),ADD_VALUES);
      	MatSetValue(*K,mapping[num_nodes+Ii],mapping[num_nodes+Ii],(*(K_tor1_holder+7)),ADD_VALUES);
      	MatSetValue(*K,mapping[num_nodes+Ii],mapping[neb_hold],(*(K_tor1_holder+8)),ADD_VALUES);
      	MatSetValue(*K,mapping[num_nodes+Ii],mapping[num_nodes+neb_hold],(*(K_tor1_holder+9)),ADD_VALUES);
      	MatSetValue(*K,mapping[num_nodes+Ii],mapping[k1],(*(K_tor1_holder+10)),ADD_VALUES);
      	MatSetValue(*K,mapping[num_nodes+Ii],mapping[num_nodes+k1],(*(K_tor1_holder+11)),ADD_VALUES);

	free(K_tor1_holder);
	K_tor1_holder=K_tor(neb_hold,Ii,k1);	    

      	MatSetValue(*K,mapping[Ii],mapping[neb_hold],(*(K_tor1_holder+12)),ADD_VALUES);
      	MatSetValue(*K,mapping[Ii],mapping[num_nodes+neb_hold],(*(K_tor1_holder+13)),ADD_VALUES);
      	MatSetValue(*K,mapping[Ii],mapping[Ii],(*(K_tor1_holder+14)),ADD_VALUES);
      	MatSetValue(*K,mapping[Ii],mapping[num_nodes+Ii],(*(K_tor1_holder+15)),ADD_VALUES);
      	MatSetValue(*K,mapping[Ii],mapping[k1],(*(K_tor1_holder+16)),ADD_VALUES);
      	MatSetValue(*K,mapping[Ii],mapping[num_nodes+k1],(*(K_tor1_holder+17)),ADD_VALUES);

      	MatSetValue(*K,mapping[num_nodes+Ii],mapping[neb_hold],(*(K_tor1_holder+18)),ADD_VALUES);
      	MatSetValue(*K,mapping[num_nodes+Ii],mapping[num_nodes+neb_hold],(*(K_tor1_holder+19)),ADD_VALUES);
      	MatSetValue(*K,mapping[num_nodes+Ii],mapping[Ii],(*(K_tor1_holder+20)),ADD_VALUES);
      	MatSetValue(*K,mapping[num_nodes+Ii],mapping[num_nodes+Ii],(*(K_tor1_holder+21)),ADD_VALUES);
      	MatSetValue(*K,mapping[num_nodes+Ii],mapping[k1],(*(K_tor1_holder+22)),ADD_VALUES);
      	MatSetValue(*K,mapping[num_nodes+Ii],mapping[num_nodes+k1],(*(K_tor1_holder+23)),ADD_VALUES);

	free(K_tor1_holder);
	
      }
    }
  }
    
  MatAssemblyBegin(*K,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*K,MAT_FINAL_ASSEMBLY);

  VecCreate(PETSC_COMM_WORLD,&F);
  PetscObjectSetName((PetscObject) F, "Load Vector");
  VecSetSizes(F,2*num_owned_nodes,PETSC_DECIDE);
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

  
  for(i=0;i<len_new;i++){
    global_surf_connectivity[2*i]=mapping[local_surf_connectivity[i]];
    global_surf_connectivity[2*i+1]=mapping[num_nodes+local_surf_connectivity[i]];
  }
  

 for(n=0;n<len_new;n++){
   x_disp[n]=phys_mesh[dimension*(local_surf_connectivity[n])] - comp_mesh[dimension*(local_surf_connectivity[n])];
   y_disp[n]=phys_mesh[dimension*(local_surf_connectivity[n])+1] - comp_mesh[dimension*(local_surf_connectivity[n])+1];
  }
 
 for(Ii=0;Ii<len_new;Ii++){
   if (local_surf_connectivity[Ii]<num_owned_nodes){
     VecSetValue(F,global_surf_connectivity[2*Ii],x_disp[Ii],INSERT_VALUES);
     VecSetValue(F,global_surf_connectivity[2*Ii+1],y_disp[Ii],INSERT_VALUES);
   }
 }
   
  VecCreate(PETSC_COMM_WORLD,&U_h);
  PetscObjectSetName((PetscObject) U_h, "lineal_torsional_smoother");
  VecSetSizes(U_h,dimension*num_owned_nodes,PETSC_DECIDE);
  VecSetFromOptions(U_h);

  PetscScalar *aU_h;
  VecGetArray(U_h,&aU_h);

  for(n=0;n<num_owned_nodes;n++){
    for(j=0;j<dimension;++j){
      aU_h[j*num_owned_nodes+n] = phys_mesh[dimension*n+j]-comp_mesh[dimension*n+j];
    }
  }
  VecRestoreArray(U_h,&aU_h);

  MatZeroRowsColumns(*K,2*len_new,global_surf_connectivity,1.0,U_h,F);

  VecAssemblyBegin(U_h);
  VecAssemblyEnd(U_h);
 
  VecAssemblyBegin(F);
  VecAssemblyEnd(F);

  MatAssemblyBegin(*K,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*K,MAT_FINAL_ASSEMBLY);

  petsc_solve_c(U_h, *K, F, options, debug_level);
  
  VecGetArray(U_h,&aU_h);
  for(n=0;n<num_owned_nodes;n++){
    for(j=0;j<dimension;++j) {
      smooth_mesh[dimension*n+j] = aU_h[j*num_owned_nodes+n]+comp_mesh[dimension*n+j];
    }
  }
  VecRestoreArray(U_h,&aU_h);
  
  VecDestroy(&F);VecDestroy(&U_h);
}
