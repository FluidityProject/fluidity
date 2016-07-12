int val(int *dim,int *n, double field[],double X[][*dim], double t) {
  int i;
  for (i=0;i<*n;i++) {
    field[i]=X[i][0]*X[i][0]-X[i][1];
  };
  return 0;
}


int val_vec(int *dim1,int *n, int* dim2,double field[][*dim2],double X[][*dim1], double t) {
  int i;
  for (i=0;i<*n;i++) {
    field[i][0]=X[i][0]*X[i][0]-X[i][1];
    field[i][1]=X[i][1]*X[i][1]-X[i][0];
  };
  return 0;
}
