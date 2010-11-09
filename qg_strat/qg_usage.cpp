/*  Copyright (C) 2006 Imperial College London and others.
    
    Please see the AUTHORS file in the main source directory for a full list
    of copyright holders.

    Prof. C Pain
    Applied Modelling and Computation Group
    Department of Earth Science and Engineering
    Imperial College London
    amcgsoftware@imperial.ac.uk
    
    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation,
    version 2.1 of the License.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
    USA
*/

#include "spud.h"
#include "qg_usage.h"

using namespace std;

using namespace Spud;

map<string, string> fl_command_line_options;

void print_version(ostream& stream){  
  stream<<"Revision: "<<__FLUIDITY_VERSION__
#ifdef DDEBUG
      <<" (debugging)"
#endif
      <<endl
      <<"Compile date: "<<__DATE__<<" "<<__TIME__<<endl
    //#include <cvsdata.h>
      <<"Adaptivity support\t\t"
#ifdef HAVE_ADAPTIVITY
      <<"yes\n"
#else
      <<"no\n"
#endif
      <<"MPI support\t\t\t"
#ifdef HAVE_MPI
      <<"yes\n"
#else
      <<"no\n"
#endif
      <<"Double precision\t\t"
#ifdef DOUBLEP
      <<"yes\n"
#else
      <<"no\n"
#endif
      <<"CGNS support\t\t\t"
#ifdef HAVE_CGNS
      <<"yes\n"
#else
      <<"no\n"
#endif
      <<"NetCDF support\t\t\t"
#ifdef HAVE_NETCDF   
      <<"yes\n"
#else
      <<"no\n"
#endif 
      <<"Signal handling support\t\t"
#ifdef NSIGNAL
      <<"no\n"
#else
      <<"yes\n"
#endif
      <<"PETSc support\t\t\t"
#ifdef HAVE_PETSC
      <<"yes\n"
#else
      <<"no\n"
#endif
      <<"ARPACK support\t\t\t"
#ifdef HAVE_LIBARPACK
      <<"yes\n"
#else
      <<"no\n"
#endif
      <<"ERA-40 forcing support\t\t"
#ifdef ENABLE_ERA40
      <<"yes\n"
#else
      <<"no\n"
#endif
      <<"Python support\t\t\t"
#ifdef HAVE_PYTHON
      <<"yes\n"
#else
      <<"no\n"
#endif
      <<"VTK support\t\t\t"
#ifdef HAVE_VTK
      <<"yes\n"
#else
      <<"no\n"
#endif
      ;
  stream.flush();
  return;
}

void print_environment(){

  char *relevant_variables[]={"FLUIDITY_INTEGER_SIZE", "FLUIDITY_REAL_SIZE",
    "PETSC_OPTIONS"};
    
  int no_relevant_variables=1;
  
  cout<<"Environment(al) variables relevant for fluidity:\n";
  for(int i=0; i<3; i++) {
    char *env=getenv(relevant_variables[i]);
    if(env!=NULL) {
      cout<<relevant_variables[i]<<" = "<<env<<endl;
      no_relevant_variables=0;
    }
  };
  if (no_relevant_variables)
    cout<<" none\n";   
  
}

void qg_usage(char *cmd){
  print_version();
  cerr<<"\n\nUsage: "<<cmd<<" [options ...] [simulation-file]\n"
      <<"\nOptions (NOTE: Long options are not available on AIX.):\n"
      <<" -h, --help\n\tHelp! Prints this message.\n"
      <<" -l, --log\n\tCreate log file for each process (useful for non-interactive testing)."
      <<" -v <level>, --verbose\n\tVerbose output to stdout, default level 0\n"
      <<" -x \n\tRun from a new xml options file\n";
  return;
}

void ParseArguments(int argc, char** argv){

#ifndef _AIX
  struct option longOptions[] = {
    {"log", 0, 0, 'l'},
    {"xml", 1, 0, 'x'},
    {"verbose", optional_argument, 0, 'v'},
    {0, 0, 0, 0}
  };
#endif
  int optionIndex = 0;
  int verbosity = 0;
  int c;

  // set opterr to nonzero to make getopt print error messages 
  opterr=1;
  
  while (true){
    c = getopt(argc, argv, "a:c:hklr:s:v::D:d:p:Vx:");
    if (c == -1) break;

    OptionError stat;
    switch (c){

    case 'h':
      fl_command_line_options["help"] = "";
      break;

    case 'l':
      fl_command_line_options["log"] = "";
      break;        

    case 'v':
      fl_command_line_options["verbose"] = (optarg == NULL) ? "1" : optarg;
      break;  
    
    case 'x':
      char xml_file[4096];
      sprintf(xml_file, "%s/%s", getenv("PWD"), optarg);
      fl_command_line_options["xml"] = xml_file;

      load_options(xml_file);
      break;
    }
  }
  
  // Help?
  if(fl_command_line_options.count("help")){
    qg_usage(argv[0]);
    exit(-1);
  }
  
  // Version?
  if(fl_command_line_options.count("version")){
    print_version();
    exit(-1);
  }
  
  // Verbose?
  {
  if(fl_command_line_options.count("verbose")==0){
    fl_command_line_options["verbose"] = "0";
  }
  
  verbosity = atoi(fl_command_line_options["verbose"].c_str());
  set_global_debug_level_fc(&verbosity);
  }
  
  // What to do with stdout?
  if(fl_command_line_options.count("log")){
  char debug_file[4096], err_file[4096];
  sprintf(debug_file, "%s/qg_strat.log", getenv("PWD"));
  sprintf(err_file, "%s/qg_strat.err", getenv("PWD"));
#ifdef HAVE_MPI
    if(MPI::Is_initialized()){
    int MyRank = MPI::COMM_WORLD.Get_rank();
    sprintf(debug_file, "%s/qg_strat.log-%d", getenv("PWD"), MyRank);
    sprintf(err_file, "%s/qg_strat.err-%d", getenv("PWD"), MyRank);
  }
#endif
  
  if(freopen(debug_file,"w",stdout)==NULL)
    perror("failed to redirect stdio for debugging"); 
    
    if(freopen(err_file,"w",stderr)==NULL)
      perror("failed to redirect stderr for debugging");
  }
  
  // Find the filename if one is specified. Assuming that final
  // argument is simulation filename if it does not correspond to a
  // known option.
  if(argc > optind + 1){
    fl_command_line_options["simulation_name"] = argv[optind + 1];
  }else if(argc == optind + 1){
    fl_command_line_options["simulation_name"] = argv[optind];
  }else if(have_option("/simulation_name")){
    OptionError stat = get_option("/simulation_name", 
                                fl_command_line_options["simulation_name"]);
    assert(stat == SPUD_NO_ERROR);
  }else{
    cerr << "ERROR: Unrecognized arguments!" << endl;
    qg_usage(argv[0]);
    exit(-1);
  }
  
  // now that we know the verbosity, and possibly redirected stdout,
  // we may print the given command line
  if(verbosity >= 2){
    cout<<"Fluidity command line:\n ";
    for(int i=0;i<argc; i++)
      cout<<argv[i]<<" ";
    cout<<endl;
    print_environment();

    // Useful for debugging options
    print_options();
  }
  
  return;
}

void PetscInit(int argc, char** argv){
#ifdef HAVE_PETSC
  int petscargc;
  char** petscargv;
  if (fl_command_line_options.count("petscopts")){
    vector<string> petscopts;
    Tokenize(fl_command_line_options["petscopts"], petscopts, ",");
    petscargc = petscopts.size() + 1;
    petscargv = new char*[petscargc];
    petscargv[0] = new char[strlen(argv[0]) + 1];
    for(size_t i = 1; i < petscargc; i++)
      petscargv[i] = new char[petscopts[i-1].size()+1];
    strncpy(petscargv[0], argv[0], strlen(argv[0]) + 1);
    for(size_t i = 1; i < petscargc; i++)
      strncpy(petscargv[i], petscopts[i-1].c_str(), petscopts[i-1].size()+1);
  }else{
    petscargc = 1;
    petscargv = new char*[1];
    petscargv[0] = new char[strlen(argv[0]) + 1];
    strncpy(petscargv[0], argv[0], strlen(argv[0]) + 1);
  }

  static char help[] = "Use --help to see the help.\n\n";
  PetscErrorCode ierr = PetscInitialize(&petscargc, &petscargv, NULL, help);
  // PetscInitializeFortran needs to be called when initialising PETSc from C, but calling it from Fortran
  // This sets all kinds of objects such as PETSC_NULL_OBJECT, PETSC_COMM_WORLD, etc., etc.
  ierr = PetscInitializeFortran();
  // CHKERRQ(ierr);
  signal(SIGSEGV, SIG_DFL);
  signal(SIGTRAP, SIG_DFL);
  
  for(size_t i = 0; i < petscargc; i++)
    delete [] petscargv[i];
  delete [] petscargv;
#endif
}

