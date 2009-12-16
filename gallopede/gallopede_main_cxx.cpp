#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <getopt.h>
#include <signal.h>

#include <map>
#include <iostream>
#include <string>

#include "fmangle.h"


extern "C" {
void gallopede_main_fc(int *inargs,int *n_args, int *u_deg);
void set_global_debug_level_fc(int *val);
}

using namespace std;

void usage(char *cmd){
  cerr<<"\n\nUsage: "<<cmd<<" [options ...]\n"
      <<"\nOptions:\n"
      <<" -h, --help\n\tHelp! Prints this message.\n"
      <<" -a, --advection\n\tTurns off advection term.\n"
      <<" -l, --linear\n\tTurns off nonlinear term.\n"
      <<" -f, --flux\n\tTurns on flux form.\n"
      <<" -p, --pressure\n\tTurns off hydrostatic pressure terms.\n"
      <<" -n, --nonhydro\n\tTurns off nonhydrostatic pressure terms.\n"
      <<" -q, --quadratic\n\tTurns on quadratic velocities.\n"
      <<" -d, --density\n\tDisables density stepping.\n"
      <<" -m, --momentum\n\tDisables momentum stepping.\n"
      <<" -v, --verbose\nSets IO dump level.\n"
      <<" -r, --rotation\n\tTurns on Coriolis terms.\n"
      <<" -g, --green\n\tTurns off non-Green-Naghdi terms.\n"
      <<" -V, --version\n\t Prints version information.\n";
  return;
}

map<string, string> fl_command_line_options;
void ParseArguments(int argc, char** argv){
  struct option longOptions[] = {
    {"help", 0, 0, 'h'},
    {"advection", 0, 0, 'a'},
    {"constant", 0, 0, 'c'},
    {"linear", 0, 0, 'l'},
    {"flux", 0, 0, 'f'},
    {"pressure", 0, 0, 'p'},
    {"quadratic", 0, 0, 'q'},
    {"density", 0, 0, 'd'},
    {"momentum", 0, 0, 'm'},
    {"upwinding", 0, 0, 'u'},
    {"verbose", 1, 0, 'v'},
    {"nonhydro", 0, 0, 'n'},
    {"nonhydro", 0, 0, 'r'},
    {"version", 0, 0, 'V'},
    {"green", 0, 0, 'g'},
    {0, 0, 0, 0}
  };
  
  int optionIndex = 0;
  int c;
  
  while (true){
    c = getopt_long(argc, argv, "haclpfqdmnuvrg::V",
		    longOptions, &optionIndex);
    if (c == -1) break;
    
    switch (c){	
    case 'h':
      fl_command_line_options["help"] = "";
      break;
      
      case 'a':
      fl_command_line_options["advection"] = "";
      break;	

    case 'c':
      fl_command_line_options["constant"] = "";
      break;	

    case 'l':
      fl_command_line_options["linear"] = "";
      break;				
      
    case 'f':
      fl_command_line_options["flux"] = "";
      break;

      case 'p':
      fl_command_line_options["pressure"] = "";
      break;

    case 'q':
      fl_command_line_options["quadratic"] = "";
      break;

 case 'd':
      fl_command_line_options["density"] = "";
      break;

 case 'm':
      fl_command_line_options["momentum"] = "";
      break;

    case 'u':
      fl_command_line_options["upwinding"] = "";
      break;

 case 'n':
      fl_command_line_options["nonhydro"] = "";
      break;

	case 'v':
      fl_command_line_options["verbose"] = (optarg == NULL) ? "3" : optarg;
      break;	

    case 'V':
      fl_command_line_options["version"] = "";
      break;

    case 'r':
      fl_command_line_options["rotation"] = "";
      break;

       case 'g':
      fl_command_line_options["green"] = "";
      break;
    }
  }
  
  return;
}

int main(int argc, char **argv){

ParseArguments(argc, argv);

if(fl_command_line_options.count("help")){
    usage(argv[0]);
    exit(-1);
  }


 {
	if(fl_command_line_options.count("verbose")==0){
	  fl_command_line_options["verbose"] = "3";
	}
	
	int level;
	level = atoi(fl_command_line_options["verbose"].c_str());
	set_global_debug_level_fc(&level);
  }


 int n_args=1;
 int inargs=0;
 int u_deg;

 u_deg=(fl_command_line_options.count("quadratic"))? 2 : 1;

if(fl_command_line_options.count("advection")){
  inargs++;
  }
 inargs=2*inargs;n_args++;
if(fl_command_line_options.count("flux")){
   inargs++;
  }
inargs=2*inargs;n_args++;
if(fl_command_line_options.count("pressure")){
    inargs++;
  }
inargs=2*inargs;n_args++;
if(fl_command_line_options.count("density")){
    inargs++;
  }
inargs=2*inargs;n_args++;
if(fl_command_line_options.count("momentum")){
    inargs++;
  }
inargs=2*inargs;n_args++;
if(fl_command_line_options.count("upwinding")){
   inargs++;
  }
inargs=2*inargs;n_args++;
if(fl_command_line_options.count("linear")){
   inargs++;
  }
inargs=2*inargs;n_args++;
if(fl_command_line_options.count("nonhydro")){
  inargs++;
  }

inargs=2*inargs;n_args++;
if(fl_command_line_options.count("constant")){
   inargs++;
  }
inargs=2*inargs;n_args++;
if(fl_command_line_options.count("rotation")){
   inargs++;
  }
inargs=2*inargs;n_args++;
if(fl_command_line_options.count("green")){
   inargs++;
  }

  gallopede_main_fc(&inargs,&n_args,&u_deg);

  exit(0);
}
