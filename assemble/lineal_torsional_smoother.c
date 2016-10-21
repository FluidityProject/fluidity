static char help[] = "2D Lineal Torsional Spring Analogy Smoother in serial and parallel.\n\n";
#include <petscksp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


void lin_tor_smoother(int dimension, int num_nodes, int num_elements, int num_surf_elements, int num_owned_nodes, int * mapping, int * connectivity, double * phys_mesh, double * smooth_mesh, double * comp_mesh, int * surf_connectivity) {
  
  Vec            F,U_h;
  Mat            K;            
  KSP            ksp;
  PC             pc;
  PetscErrorCode ierr;
  PetscInt       i,n,m,j,its,Ii,len_new=0,local_surf_connectivity[dimension*num_surf_elements], global_surf_connectivity[dimension*num_surf_elements],
                 num_nodes_col_x[num_nodes],num_nodes_col_y[num_nodes];
  PetscScalar    length,x_disp[num_surf_elements],y_disp[num_surf_elements],z_disp[num_surf_elements],smoothed_x[num_nodes],smoothed_y[num_nodes],smoothed_z[num_nodes];


  int conn_mat[num_elements][3];
  for(m=0;m<num_elements;m++){
    conn_mat[m][0]=*(connectivity+(3*m));
    conn_mat[m][1]=*(connectivity+(3*m)+1);
    conn_mat[m][2]=*(connectivity+(3*m)+2);
  }

  double lij(int dimension, int nodei, int nodej){
    nodei--;nodej--;
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

  double alpha_ang(int nodei, int nodej){
    nodei--;nodej--;
    double alpha;
    double xi = comp_mesh[dimension*nodei]; double yi = comp_mesh[dimension*nodei+1];
    double xj = comp_mesh[dimension*nodej]; double yj = comp_mesh[dimension*nodej+1];
    double x_edge = xj - xi; double y_edge = yj - yi;
    alpha = atan2(y_edge,x_edge);
    return alpha;
  }
    
  /*Lineal Stiffness*/
  double * K_lin(int nodei, int nodej){
    static double K_lin_temp[16];
    double inv_length = 1.0/lij(2,nodei,nodej);
    double alpha = alpha_ang(nodei,nodej);
    
    K_lin_temp[0]=inv_length*cos(alpha)*cos(alpha);K_lin_temp[1]=inv_length*sin(alpha)*cos(alpha);K_lin_temp[2]=-1.0*inv_length*cos(alpha)*cos(alpha);K_lin_temp[3]=-1.0*inv_length*cos(alpha)*sin(alpha);
    K_lin_temp[4]=inv_length*sin(alpha)*cos(alpha);K_lin_temp[5]=inv_length*sin(alpha)*sin(alpha);K_lin_temp[6]=-1.0*inv_length*sin(alpha)*cos(alpha);K_lin_temp[7]=-1.0*inv_length*sin(alpha)*sin(alpha);
    K_lin_temp[8]=-1.0*inv_length*cos(alpha)*cos(alpha);K_lin_temp[9]=-1.0*inv_length*sin(alpha)*cos(alpha);K_lin_temp[10]=inv_length*cos(alpha)*cos(alpha);K_lin_temp[11]=inv_length*sin(alpha)*cos(alpha);
    K_lin_temp[12]=-1.0*inv_length*sin(alpha)*cos(alpha);K_lin_temp[13]=-1.0*inv_length*sin(alpha)*sin(alpha);K_lin_temp[14]=inv_length*sin(alpha)*cos(alpha);K_lin_temp[15]=inv_length*sin(alpha)*sin(alpha);

    return K_lin_temp;
  }

  double area(int nodei, int nodej, int nodek){
    double s = (lij(2,nodei,nodej) + lij(2,nodei,nodek) + lij(2,nodej,nodek))/2.0;
    double area = sqrt(s*(s-lij(2,nodei,nodej))*(s-lij(2,nodei,nodek))*(s-lij(2,nodej,nodek)));
    return area;
  }

  double tors(int nodei, int nodej, int nodek){
    double tors = (pow(lij(2,nodei,nodej),2) * pow(lij(2,nodei,nodek),2))/(4*pow(area(nodei,nodej,nodek),2));
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
    nodei--;nodej--;
    xi = comp_mesh[dimension*nodei];
    xj = comp_mesh[dimension*nodej];
    x_edge = xj - xi;
    return x_edge;
  }

  double y_edge(int nodei,int nodej){
    double yi,yj,y_edge;
    nodei--;nodej--;
    yi = comp_mesh[dimension*nodei+1];
    yj = comp_mesh[dimension*nodej+1];
    y_edge = yj - yi;
    return y_edge;
  }

  double a_rot(int nodei, int nodej){
    double a_rot = x_edge(nodei,nodej)/(pow(lij(2,nodei,nodej),2));
    return a_rot;
  }

  double b_rot(int nodei, int nodej){
    double b_rot = y_edge(nodei,nodej)/(pow(lij(2,nodei,nodej),2));
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
    int i,j,k,r1=6,c1=3,r2,c2=3;
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
    
    r1=6;c1=3;r2=3;c2=6;
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

  int *neb_tri(int nodei, int nodej,int * connectivity){
    int n,i=0;
    static int k[2];

    while(i<3){
    for(n=0;n<num_elements;n++){
      if(nodei==*(connectivity+(3*n))){
	if(nodej==*(connectivity+(3*n)+1)){
	  k[i]=*(connectivity+(3*n)+2);
	  i++;
	  continue;
	   }
      }
      if(nodei==*(connectivity+(3*n))){
	if(nodej==*(connectivity+(3*n)+2)){
	  k[i]=*(connectivity+(3*n)+1);
	  i++;
	  continue;
	}
      }
      if(nodei==*(connectivity+(3*n)+1)){
	if(nodej==*(connectivity+(3*n))){
	  k[i]=*(connectivity+(3*n)+2);
	  i++;
	  continue;
	}
      }
      if(nodei==*(connectivity+(3*n)+1)){
	if(nodej==*(connectivity+(3*n)+2)){
	  k[i]=*(connectivity+(3*n));
	  i++;
	  continue;
	}
      }
      if(nodei==*(connectivity+(3*n)+2)){
	if(nodej==*(connectivity+(3*n))){
	  k[i]=*(connectivity+(3*n)+1);
	  i++;
	  continue;
	}
      }
      if(nodei==*(connectivity+(3*n)+2)){
	if(nodej==*(connectivity+(3*n)+1)){
	  k[i]=*(connectivity+(3*n));
	  i++;
	  continue;
	}
      }
    }
    }
    if(k[0]==k[1]){
      k[1]=0;}
    return k;
  }

  /*Neighbour Matrix*/
  int neb_mat[num_owned_nodes][8];
  for(Ii=0;Ii<num_owned_nodes;Ii++){
    int neb_a_count=0; int neb_b_count=0; int neb_c_count=0;
    int neb_a[8]={};int neb_b[8]={};int neb_c[8]={};int neb_tot[24]={};int neb_fin[8]={};
    
    for(m=0;m<num_elements;++m){
     
      if (Ii+1 == conn_mat[m][0]){
	neb_a[neb_a_count]=conn_mat[m][1];
	neb_a[neb_a_count+1]=conn_mat[m][2];
	neb_a_count+=2;}
      
      if (Ii+1 == conn_mat[m][1]){
	neb_b[neb_b_count]=conn_mat[m][0];
	neb_b[neb_b_count+1]=conn_mat[m][2];
	neb_b_count+=2;}
      
      if (Ii+1 == conn_mat[m][2]){
	  neb_c[neb_c_count]=conn_mat[m][0];
	  neb_c[neb_c_count+1]=conn_mat[m][1];
	  neb_c_count+=2;}
    }
    
    for(n=0;n<8;n++){
      neb_tot[n]=neb_a[n];
      neb_tot[8+n]=neb_b[n];
      neb_tot[16+n]=neb_c[n];
    }

   //Removing duplicates 
   for(m=0;m<24;++m){
     for(n=m+1;n<24;++n){
       if(neb_tot[m]==neb_tot[n]){
    	  neb_tot[n]=0;
       }
     }
   }

   //sorting
   for(i=1; i<=24-1; i++){
     for(j=1; j<=24-i; j++){
       if(neb_tot[j-1] <= neb_tot[j]){
	 int t = neb_tot[j-1];
	 neb_tot[j-1]=neb_tot[j];
	 neb_tot[j]=t;
       }
     }
   }

   //Only interested in first 8 ints
   for(m=0;m<8;m++){
     neb_mat[Ii][m]=neb_tot[m];}
  }

  MatCreate(PETSC_COMM_WORLD,&K);
  PetscObjectSetName((PetscObject) K, "Stiffness Matrix");
  MatSetSizes(K,2*num_owned_nodes,2*num_owned_nodes,PETSC_DETERMINE,PETSC_DETERMINE);
  MatSetUp(K);


  for(Ii=0;Ii<num_owned_nodes;Ii++){
       for(n=0;n<8;n++){
     	  int neb_hold = neb_mat[Ii][n];
	  if(neb_hold != 0){
	    
	    MatSetValue(K,2*mapping[Ii],2*mapping[Ii],*(K_lin(Ii+1,neb_hold)+0),ADD_VALUES);
	    MatSetValue(K,2*mapping[Ii],2*mapping[Ii]+1,*(K_lin(Ii+1,neb_hold)+1),ADD_VALUES);
	    MatSetValue(K,2*mapping[Ii],2*mapping[neb_hold-1],*(K_lin(Ii+1,neb_hold)+2),ADD_VALUES);
	    MatSetValue(K,2*mapping[Ii],2*mapping[neb_hold-1]+1,*(K_lin(Ii+1,neb_hold)+3),ADD_VALUES);
	    
	    MatSetValue(K,2*mapping[Ii]+1,2*mapping[Ii],*(K_lin(Ii+1,neb_hold)+4),ADD_VALUES);
	    MatSetValue(K,2*mapping[Ii]+1,2*mapping[Ii]+1,*(K_lin(Ii+1,neb_hold)+5),ADD_VALUES);
	    MatSetValue(K,2*mapping[Ii]+1,2*mapping[neb_hold-1],*(K_lin(Ii+1,neb_hold)+6),ADD_VALUES);
	    MatSetValue(K,2*mapping[Ii]+1,2*mapping[neb_hold-1]+1,*(K_lin(Ii+1,neb_hold)+7),ADD_VALUES);

	    MatSetValue(K,2*mapping[neb_hold-1],2*mapping[Ii],*(K_lin(Ii+1,neb_hold)+8),ADD_VALUES);
	    MatSetValue(K,2*mapping[neb_hold-1],2*mapping[Ii]+1,*(K_lin(Ii+1,neb_hold)+9),ADD_VALUES);
	    MatSetValue(K,2*mapping[neb_hold-1],2*mapping[neb_hold-1],*(K_lin(Ii+1,neb_hold)+10),ADD_VALUES);
	    MatSetValue(K,2*mapping[neb_hold-1],2*mapping[neb_hold-1]+1,*(K_lin(Ii+1,neb_hold)+11),ADD_VALUES);
	    
	    MatSetValue(K,2*mapping[neb_hold-1]+1,2*mapping[Ii],*(K_lin(Ii+1,neb_hold)+12),ADD_VALUES);
	    MatSetValue(K,2*mapping[neb_hold-1]+1,2*mapping[Ii]+1,*(K_lin(Ii+1,neb_hold)+13),ADD_VALUES);
	    MatSetValue(K,2*mapping[neb_hold-1]+1,2*mapping[neb_hold-1],*(K_lin(Ii+1,neb_hold)+14),ADD_VALUES);
	    MatSetValue(K,2*mapping[neb_hold-1]+1,2*mapping[neb_hold-1]+1,*(K_lin(Ii+1,neb_hold)+15),ADD_VALUES);

	    
	    int k1=*(neb_tri(Ii+1,neb_hold,connectivity)+0);int k2=*(neb_tri(Ii+1,neb_hold,connectivity)+1);
	    double *K_tor1_holder=K_tor(Ii+1,neb_hold,k1);	    

	    MatSetValue(K,2*mapping[Ii],2*mapping[Ii],*(K_tor1_holder+0),ADD_VALUES);
	    MatSetValue(K,2*mapping[Ii],2*mapping[Ii]+1,*(K_tor1_holder+1),ADD_VALUES);
	    MatSetValue(K,2*mapping[Ii],2*mapping[neb_hold-1],*(K_tor1_holder+2),ADD_VALUES);
	    MatSetValue(K,2*mapping[Ii],2*mapping[neb_hold-1]+1,*(K_tor1_holder+3),ADD_VALUES);
	    MatSetValue(K,2*mapping[Ii],2*mapping[k1-1],*(K_tor1_holder+4),ADD_VALUES);
	    MatSetValue(K,2*mapping[Ii],2*mapping[k1-1]+1,*(K_tor1_holder+5),ADD_VALUES);

	    MatSetValue(K,2*mapping[Ii]+1,2*mapping[Ii],*(K_tor1_holder+6),ADD_VALUES);
	    MatSetValue(K,2*mapping[Ii]+1,2*mapping[Ii]+1,*(K_tor1_holder+7),ADD_VALUES);
	    MatSetValue(K,2*mapping[Ii]+1,2*mapping[neb_hold-1],*(K_tor1_holder+8),ADD_VALUES);
	    MatSetValue(K,2*mapping[Ii]+1,2*mapping[neb_hold-1]+1,*(K_tor1_holder+9),ADD_VALUES);
	    MatSetValue(K,2*mapping[Ii]+1,2*mapping[k1-1],*(K_tor1_holder+10),ADD_VALUES);
	    MatSetValue(K,2*mapping[Ii]+1,2*mapping[k1-1]+1,*(K_tor1_holder+11),ADD_VALUES);

	    MatSetValue(K,2*mapping[neb_hold-1],2*mapping[Ii],*(K_tor1_holder+12),ADD_VALUES);
	    MatSetValue(K,2*mapping[neb_hold-1],2*mapping[Ii]+1,*(K_tor1_holder+13),ADD_VALUES);
	    MatSetValue(K,2*mapping[neb_hold-1],2*mapping[neb_hold-1],*(K_tor1_holder+14),ADD_VALUES);
	    MatSetValue(K,2*mapping[neb_hold-1],2*mapping[neb_hold-1]+1,*(K_tor1_holder+15),ADD_VALUES);
	    MatSetValue(K,2*mapping[neb_hold-1],2*mapping[k1-1],*(K_tor1_holder+16),ADD_VALUES);
	    MatSetValue(K,2*mapping[neb_hold-1],2*mapping[k1-1]+1,*(K_tor1_holder+17),ADD_VALUES);

	    MatSetValue(K,2*mapping[neb_hold-1]+1,2*mapping[Ii],*(K_tor1_holder+18),ADD_VALUES);
	    MatSetValue(K,2*mapping[neb_hold-1]+1,2*mapping[Ii]+1,*(K_tor1_holder+19),ADD_VALUES);
            MatSetValue(K,2*mapping[neb_hold-1]+1,2*mapping[neb_hold-1],*(K_tor1_holder+20),ADD_VALUES);
	    MatSetValue(K,2*mapping[neb_hold-1]+1,2*mapping[neb_hold-1]+1,*(K_tor1_holder+21),ADD_VALUES);
	    MatSetValue(K,2*mapping[neb_hold-1]+1,2*mapping[k1-1],*(K_tor1_holder+22),ADD_VALUES);
	    MatSetValue(K,2*mapping[neb_hold-1]+1,2*mapping[k1-1]+1,*(K_tor1_holder+23),ADD_VALUES);

	    free(K_tor1_holder);
	    
	    if(k2 != 0){
	      double *K_tor2_holder=K_tor(Ii+1,neb_hold,k2);

	      MatSetValue(K,2*mapping[Ii],2*mapping[Ii],*(K_tor2_holder+0),ADD_VALUES);
	      MatSetValue(K,2*mapping[Ii],2*mapping[Ii]+1,*(K_tor2_holder+1),ADD_VALUES);
	      MatSetValue(K,2*mapping[Ii],2*mapping[neb_hold-1],*(K_tor2_holder+2),ADD_VALUES);
	      MatSetValue(K,2*mapping[Ii],2*mapping[neb_hold-1]+1,*(K_tor2_holder+3),ADD_VALUES);
	      MatSetValue(K,2*mapping[Ii],2*mapping[k2-1],*(K_tor2_holder+4),ADD_VALUES);
	      MatSetValue(K,2*mapping[Ii],2*mapping[k2-1]+1,*(K_tor2_holder+5),ADD_VALUES);

	      MatSetValue(K,2*mapping[Ii]+1,2*mapping[Ii],*(K_tor2_holder+6),ADD_VALUES);
	      MatSetValue(K,2*mapping[Ii]+1,2*mapping[Ii]+1,*(K_tor2_holder+7),ADD_VALUES);
	      MatSetValue(K,2*mapping[Ii]+1,2*mapping[neb_hold-1],*(K_tor2_holder+8),ADD_VALUES);
	      MatSetValue(K,2*mapping[Ii]+1,2*mapping[neb_hold-1]+1,*(K_tor2_holder+9),ADD_VALUES);
	      MatSetValue(K,2*mapping[Ii]+1,2*mapping[k2-1],*(K_tor2_holder+10),ADD_VALUES);
	      MatSetValue(K,2*mapping[Ii]+1,2*mapping[k2-1]+1,*(K_tor2_holder+11),ADD_VALUES);

	      MatSetValue(K,2*mapping[neb_hold-1],2*mapping[Ii],*(K_tor2_holder+12),ADD_VALUES);
	      MatSetValue(K,2*mapping[neb_hold-1],2*mapping[Ii]+1,*(K_tor2_holder+13),ADD_VALUES);
	      MatSetValue(K,2*mapping[neb_hold-1],2*mapping[neb_hold-1],*(K_tor2_holder+14),ADD_VALUES);
	      MatSetValue(K,2*mapping[neb_hold-1],2*mapping[neb_hold-1]+1,*(K_tor2_holder+15),ADD_VALUES);
	      MatSetValue(K,2*mapping[neb_hold-1],2*mapping[k2-1],*(K_tor2_holder+16),ADD_VALUES);
	      MatSetValue(K,2*mapping[neb_hold-1],2*mapping[k2-1]+1,*(K_tor2_holder+17),ADD_VALUES);

	      MatSetValue(K,2*mapping[neb_hold-1]+1,2*mapping[Ii],*(K_tor2_holder+18),ADD_VALUES);
	      MatSetValue(K,2*mapping[neb_hold-1]+1,2*mapping[Ii]+1,*(K_tor2_holder+19),ADD_VALUES);
	      MatSetValue(K,2*mapping[neb_hold-1]+1,2*mapping[neb_hold-1],*(K_tor2_holder+20),ADD_VALUES);
	      MatSetValue(K,2*mapping[neb_hold-1]+1,2*mapping[neb_hold-1]+1,*(K_tor2_holder+21),ADD_VALUES);
	      MatSetValue(K,2*mapping[neb_hold-1]+1,2*mapping[k2-1],*(K_tor2_holder+22),ADD_VALUES);
	      MatSetValue(K,2*mapping[neb_hold-1]+1,2*mapping[k2-1]+1,*(K_tor2_holder+23),ADD_VALUES);
	      
	      free(K_tor2_holder);
	    }
	    
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
    global_surf_connectivity[2*i]=2*mapping[local_surf_connectivity[i]];
    global_surf_connectivity[2*i+1]=2*mapping[local_surf_connectivity[i]]+1;
  }
  

  MatZeroRows(K,2*len_new,global_surf_connectivity,1.0,NULL,NULL);

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
