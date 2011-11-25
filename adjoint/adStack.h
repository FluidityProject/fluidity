
#ifndef ADSTACK_LOADED
#define ADSTACK_LOADED 1

extern void pushcharacterarray(char *x, int n) ;
extern void popcharacterarray(char *x, int n) ;
extern void lookcharacterarray(char *x, int n) ;

extern void pushbooleanarray(char *x, int n) ;
extern void popbooleanarray(char *x, int n) ;
extern void lookbooleanarray(char *x, int n) ;

extern void pushinteger4array(int *x, int n) ;
extern void popinteger4array(int *x, int n) ;
extern void lookinteger4array(int *x, int n) ;

extern void pushinteger8array(long int *x, int n) ;
extern void popinteger8array(long int *x, int n) ;
extern void lookinteger8array(long int *x, int n) ;

extern void pushinteger16array(long long int *x, int n) ;
extern void popinteger16array(long long int *x, int n) ;
extern void lookinteger16array(long long int *x, int n) ;

extern void pushreal4array(float *x, int n) ;
extern void popreal4array(float *x, int n) ;
extern void lookreal4array(float *x, int n) ;

extern void pushreal8array(double *x, int n) ;
extern void popreal8array(double *x, int n) ;
extern void lookreal8array(double *x, int n) ;

extern void pushreal16array(char *x, int n) ;
extern void popreal16array(char *x, int n) ;
extern void lookreal16array(char *x, int n) ;

extern void pushreal32array(char *x, int n) ;
extern void popreal32array(char *x, int n) ;
extern void lookreal32array(char *x, int n) ;

extern void pushcomplex4array(char *x, int n) ;
extern void popcomplex4array(char *x, int n) ;
extern void lookcomplex4array(char *x, int n) ;

extern void pushcomplex8array(char *x, int n) ;
extern void popcomplex8array(char *x, int n) ;
extern void lookcomplex8array(char *x, int n) ;

extern void pushcomplex16array(char *x, int n) ;
extern void popcomplex16array(char *x, int n) ;
extern void lookcomplex16array(char *x, int n) ;

extern void pushcomplex32array(char *x, int n) ;
extern void popcomplex32array(char *x, int n) ;
extern void lookcomplex32array(char *x, int n) ;

extern void pushNarray(char *x, unsigned int nbChars) ;
extern void popNarray(char *x, unsigned int nbChars) ;
extern void lookNarray(char *x, unsigned int nbChars) ;

extern void resetADLookStack() ;

extern void printbigbytes(long int nbblocks, long int blocksz, long int nbunits) ;

extern void printctraffic_() ;

extern void printtopplace_() ;

extern void printstackmax_() ;

extern void printlookingplace_() ;

extern void showrecentcstack_() ;

#endif
