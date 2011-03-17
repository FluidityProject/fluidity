/* A cpp main for calling MultiphasePrototype */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>

#define F77_FUNC(name,NAME) name ## _
#define F77_FUNC_(name,NAME) name ## _

extern "C" {
#define multiphase_prototype_fc F77_FUNC(multiphase_prototype, MULTIPHASE_PROTOTYPE)
  void multiphase_prototype_fc();
#ifdef USING_GFORTRAN
  /* gfortran hack to ensure 4-byte record marker for unformatted files */
  void _gfortran_set_record_marker(int);
#endif
}


int main(){
  
#ifdef USING_GFORTRAN
  /* gfortran hack to ensure 4-byte record marker for unformatted files */
  _gfortran_set_record_marker(4);
#endif
     
  // Start simulator
  multiphase_prototype_fc();    

  exit(0);
}









