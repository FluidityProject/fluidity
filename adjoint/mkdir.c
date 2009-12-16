#include <sys/stat.h>
#include <sys/types.h>
#include <stdlib.h>
#include <string.h>

#include "fmangle.h"

void flmkdir_fc(char *path, int *len, int *ierr){
  char *dir;
  dir = (char *)malloc(*len + 1);
  strncpy(dir, path, *len);
  dir[*len]='\0';
  
  *ierr = mkdir(dir, 0700);
  
  return;
}

