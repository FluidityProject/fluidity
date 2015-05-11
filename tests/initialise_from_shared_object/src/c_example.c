int val(int *dim,int *n, double field[],double X[][*dim], double t) {
  int i;
  for (i=0;i<*n;i++) {
    field[i]=X[i][0]*X[i][0]-X[i][1];
  };
  return 0;
}
