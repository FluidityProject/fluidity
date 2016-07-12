extern "C" int val(int*,int*,double *,void *,double);
extern "C" int val_vec(int*,int*,int*,void *,void *,double);

int val(int *dim,int *n, double field[],void *X2d, double t) {
  // double (*X)[*dim]= (double (*)[*dim])(X2d);
  for (int i=0;i<*n;i++) {
    field[i]=((double (*)[*dim])(X2d))[i][0]*((double (*)[*dim])(X2d))[i][0]-((double (*)[*dim])(X2d))[i][1];
  };
  return 0;
}



int val_vec(int *dim1,int *n,int *dim2, void *field,void *X2d, double t) {
  // double (*X)[*dim]= (double (*)[*dim])(X2d);
  for (int i=0;i<*n;i++) {
    ((double (*)[*dim2])field)[i][0]=((double (*)[*dim1])(X2d))[i][0]*((double (*)[*dim1])(X2d))[i][0]-((double (*)[*dim1])(X2d))[i][1];
    ((double (*)[*dim2])field)[i][1]=((double (*)[*dim1])(X2d))[i][1]*((double (*)[*dim1])(X2d))[i][1]-((double (*)[*dim1])(X2d))[i][0];
  };
  return 0;
}
