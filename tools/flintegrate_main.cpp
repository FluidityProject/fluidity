#include "confdefs.h"
#include <iostream>
#include <stdlib.h>
#include <string.h>

#define partial_integrals_fc F77_FUNC_(partial_integrals, PARTIAL_INTEGRALS)
extern "C" {
  void partial_integrals_fc(const char *simulationname, int *simulationname_len,
                            const char *dumpnumber, int *dumpnumber_len,
                            const char *vol_tracer_name, int *vol_tracer_name_len,
                            const char *int_tracer_name, int *int_tracer_name_len);
}

using namespace std;

int main(int argc, char *argv[]){
  // Reall should have a test and usage function, but it wasn't in the
  // original code...
  if(argc!=5){
    cerr<<"ERROR: wrong number of arguments\n";
    exit(-1);
  }

  int len[4];
  for(size_t i=0; i<4; i++){
    len[i] = strlen(argv[i+1]);
  }
  
  partial_integrals_fc(argv[1], len,
                       argv[2], len+1,
                       argv[3], len+2,
                       argv[4], len+3);

  return 0;
}
