#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <getopt.h>
#include <signal.h>

#include <map>
#include <iostream>
#include <string>

#include "fmangle.h"

using namespace std;

int main(int argc, char **argv){

  int level = 3;

  set_global_debug_level_fc(&level);
  gn_main_fc();

  exit(0);
}
