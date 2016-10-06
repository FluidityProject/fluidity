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

  
  //Euclidean length function for [1,2,3] space
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
      double xi = comp_mesh[dimension*nodei-2]; double yi = comp_mesh[dimension*nodei-1];
      double xj = comp_mesh[dimension*nodej-2]; double yj = comp_mesh[dimension*nodej-1];
      length = sqrt(pow(xi-xj,2)+pow(yi-yj,2));
    }
    else {
      double xi = comp_mesh[dimension*nodei-3]; double yi = comp_mesh[dimension*nodei-2]; double zi = comp_mesh[dimension*nodei-1];
      double xj = comp_mesh[dimension*nodej-3]; double yj = comp_mesh[dimension*nodej-2]; double zj = comp_mesh[dimension*nodej-1];
      length = sqrt(pow(xi-xj,2)+pow(yi-yj,2)+pow(zi-zj,2));
    }
    
    return length;
  }
  

  double * K_lin(int nodei, int nodej){
    static double K_lin_temp[16];
    double xi = comp_mesh[dimension*nodei-2]; double yi = comp_mesh[dimension*nodei-1];
    double xj = comp_mesh[dimension*nodej-2]; double yj = comp_mesh[dimension*nodej-1];
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
    conn_mat[m][0]=mapping[*(connectivity+(3*m))];
    conn_mat[m][1]=mapping[*(connectivity+(3*m)+1)];
    conn_mat[m][2]=mapping[*(connectivity+(3*m)+2)];
    PetscPrintf(PETSC_COMM_SELF,"conn_mat[%d]=[%d,%d,%d]\n",m,conn_mat[m][0],conn_mat[m][1],conn_mat[m][2]);
  }
  
  int * neb(int nodei){
    int neb_a_count=0; int neb_b_count=0; int neb_c_count=0;
    int neb_a[8]={};int neb_b[8]={};int neb_c[8]={};int neb_tot[24]={};

    int* tot_fin = malloc(8*sizeof(int));
    int* tot = malloc(24*sizeof(int));
    
    for(int m=0;m<num_elements;++m){
      if (nodei == conn_mat[m][0]){
	neb_a[neb_a_count]=conn_mat[m][1];
	neb_a[neb_a_count+1]=conn_mat[m][2];
	neb_a_count+=2;
      }
      if (nodei == conn_mat[m][1]){
	neb_b[neb_b_count]=conn_mat[m][0];
	neb_b[neb_b_count+1]=conn_mat[m][2];
	neb_b_count+=2;
      }
      if (nodei == conn_mat[m][2]){
	  neb_c[neb_c_count]=conn_mat[m][0];
	  neb_c[neb_c_count]=conn_mat[m][1];
	  neb_c_count+=2;
	}
    }
	  
      memcpy(tot,neb_a,8*sizeof(int));
      memcpy(tot+8,neb_b,8*sizeof(int));
      memcpy(tot+16,neb_c,8*sizeof(int));

     //removing duplicates 
     for(int m=0;m<24;++m){
   	for(int n=m+1;n<24;++n){
    	  if(tot[m]==tot[n]){
    	    tot[n]=0;
    	  }
    	} 
     }

       //sorting
        for(int i=1; i<=24-1; i++){                               
	 for(int j=1; j<=24-i; j++){                             
	   if(tot[j-1] <= tot[j]){
	       int t = tot[j-1];
	       tot[j-1]=tot[j];
	       tot[j]=t;
             }
	 }
       }
	
	memcpy(tot_fin,tot,8*sizeof(int));
	free(tot);
   
    return tot_fin;

  }

  
  MatCreate(PETSC_COMM_WORLD,&K);
  PetscObjectSetName((PetscObject) K, "Stiffness Matrix");
  MatSetSizes(K,2*num_owned_nodes,2*num_owned_nodes,PETSC_DETERMINE,PETSC_DETERMINE);
  MatSetUp(K);

  MatGetOwnershipRange(K,&Istart,&Iend);
  PetscPrintf(PETSC_COMM_SELF,"num_owned_nodes=%d\n",num_owned_nodes);

  
  for(int Ii=0;Ii<num_nodes;Ii++){
    for(int n=0;n<8;n++){
      if(Ii<=num_owned_nodes){ 
	int neb_hold = *(neb(Ii)+n);
	/*    if(neb_hold != 0){
	
	MatSetValue(K,2*(mapping[Ii]+1)-2,2*(mapping[Ii]+1)-2,*(K_lin(Ii,neb_hold)+0),ADD_VALUES);
	MatSetValue(K,2*(mapping[Ii]+1)-2,2*(mapping[Ii]+1)-1,*(K_lin(Ii,neb_hold)+1),ADD_VALUES);
	MatSetValue(K,2*(mapping[Ii]+1)-2,2*(mapping[neb_hold]+1)-2,*(K_lin(Ii,neb_hold)+2),ADD_VALUES);
	MatSetValue(K,2*(mapping[Ii]+1)-2,2*(mapping[neb_hold]+1)-1,*(K_lin(Ii,neb_hold)+3),ADD_VALUES);

	MatSetValue(K,2*(mapping[Ii]+1)-1,2*(mapping[Ii]+1)-2,*(K_lin(Ii,neb_hold)+4),ADD_VALUES);
	MatSetValue(K,2*(mapping[Ii]+1)-1,2*(mapping[Ii]+1)-1,*(K_lin(Ii,neb_hold)+5),ADD_VALUES);
	MatSetValue(K,2*(mapping[Ii]+1)-1,2*(mapping[neb_hold]+1)-2,*(K_lin(Ii,neb_hold)+6),ADD_VALUES);
	MatSetValue(K,2*(mapping[Ii]+1)-1,2*(mapping[neb_hold]+1)-1,*(K_lin(Ii,neb_hold)+7),ADD_VALUES);

	MatSetValue(K,2*(mapping[neb_hold]+1)-2,2*(mapping[Ii]+1)-2,*(K_lin(Ii,neb_hold)+8),ADD_VALUES);
	MatSetValue(K,2*(mapping[neb_hold]+1)-2,2*(mapping[Ii]+1)-1,*(K_lin(Ii,neb_hold)+9),ADD_VALUES);
	MatSetValue(K,2*(mapping[neb_hold]+1)-2,2*(mapping[neb_hold]+1)-2,*(K_lin(Ii,neb_hold)+10),ADD_VALUES);
	MatSetValue(K,2*(mapping[neb_hold]+1)-2,2*(mapping[neb_hold]+1)-1,*(K_lin(Ii,neb_hold)+11),ADD_VALUES);

	MatSetValue(K,2*(mapping[neb_hold]+1)-1,2*(mapping[Ii]+1)-2,*(K_lin(Ii,neb_hold)+12),ADD_VALUES);
	MatSetValue(K,2*(mapping[neb_hold]+1)-1,2*(mapping[Ii]+1)-1,*(K_lin(Ii,neb_hold)+13),ADD_VALUES);
	MatSetValue(K,2*(mapping[neb_hold]+1)-1,2*(mapping[neb_hold]+1)-2,*(K_lin(Ii,neb_hold)+14),ADD_VALUES);
	MatSetValue(K,2*(mapping[neb_hold]+1)-1,2*(mapping[neb_hold]+1)-1,*(K_lin(Ii,neb_hold)+15),ADD_VALUES);
      
	
	MatSetValue(K,2*(mapping[Ii-1]+1)-2,2*(mapping[Ii-1]+1)-2,*(K_lin(Ii,neb_hold)+0),ADD_VALUES);
	MatSetValue(K,2*(mapping[Ii-1]+1)-2,2*(mapping[Ii-1]+1)-1,*(K_lin(Ii,neb_hold)+1),ADD_VALUES);
	MatSetValue(K,2*(mapping[Ii-1]+1)-2,2*(mapping[neb_hold-1]+1)-2,*(K_lin(Ii,neb_hold)+2),ADD_VALUES);
	MatSetValue(K,2*(mapping[Ii-1]+1)-2,2*(mapping[neb_hold-1]+1)-1,*(K_lin(Ii,neb_hold)+3),ADD_VALUES);

	MatSetValue(K,2*(mapping[Ii-1]+1)-1,2*(mapping[Ii-1]+1)-2,*(K_lin(Ii,neb_hold)+4),ADD_VALUES);
	MatSetValue(K,2*(mapping[Ii-1]+1)-1,2*(mapping[Ii-1]+1)-1,*(K_lin(Ii,neb_hold)+5),ADD_VALUES);
	MatSetValue(K,2*(mapping[Ii-1]+1)-1,2*(mapping[neb_hold-1]+1)-2,*(K_lin(Ii,neb_hold)+6),ADD_VALUES);
	MatSetValue(K,2*(mapping[Ii-1]+1)-1,2*(mapping[neb_hold-1]+1)-1,*(K_lin(Ii,neb_hold)+7),ADD_VALUES);

	MatSetValue(K,2*(mapping[neb_hold-1]+1)-2,2*(mapping[Ii-1]+1)-2,*(K_lin(Ii,neb_hold)+8),ADD_VALUES);
	MatSetValue(K,2*(mapping[neb_hold-1]+1)-2,2*(mapping[Ii-1]+1)-1,*(K_lin(Ii,neb_hold)+9),ADD_VALUES);
	MatSetValue(K,2*(mapping[neb_hold-1]+1)-2,2*(mapping[neb_hold-1]+1)-2,*(K_lin(Ii,neb_hold)+10),ADD_VALUES);
	MatSetValue(K,2*(mapping[neb_hold-1]+1)-2,2*(mapping[neb_hold-1]+1)-1,*(K_lin(Ii,neb_hold)+11),ADD_VALUES);

	MatSetValue(K,2*(mapping[neb_hold-1]+1)-1,2*(mapping[Ii-1]+1)-2,*(K_lin(Ii,neb_hold)+12),ADD_VALUES);
	MatSetValue(K,2*(mapping[neb_hold-1]+1)-1,2*(mapping[Ii-1]+1)-1,*(K_lin(Ii,neb_hold)+13),ADD_VALUES);
	MatSetValue(K,2*(mapping[neb_hold-1]+1)-1,2*(mapping[neb_hold-1]+1)-2,*(K_lin(Ii,neb_hold)+14),ADD_VALUES);
	MatSetValue(K,2*(mapping[neb_hold-1]+1)-1,2*(mapping[neb_hold-1]+1)-1,*(K_lin(Ii,neb_hold)+15),ADD_VALUES);
		
	}*/
    }
    }
  }
   
 
  
  // PetscPrintf(PETSC_COMM_WORLD,"exited");
  
  VecCreate(PETSC_COMM_WORLD,&F);
  PetscObjectSetName((PetscObject) F, "Load Vector");
  VecSetSizes(F,2*num_owned_nodes,PETSC_DECIDE);
  VecSetFromOptions(F);
  
  MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);

  VecAssemblyBegin(F);
  VecAssemblyEnd(F);

// MatView(K,PETSC_VIEWER_DEFAULT);
 
  VecDestroy(&F);
  MatDestroy(&K);
 }
