/*  Copyright (C) 2006 Imperial College London and others.
    
    Please see the AUTHORS file in the main source directory for a full list
    of copyright holders.

    Prof. C Pain
    Applied Modelling and Computation Group
    Department of Earth Science and Engineering
    Imperial College London

    C.Pain@Imperial.ac.uk
    
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

#include "Usage.h"
#include "spud.h"

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
      <<"2D adaptivity support\t\t"
#ifdef HAVE_MBA_2D
      <<"yes\n"
#else
      <<"no\n"
#endif
      <<"3D MBA support\t\t\t"
#ifdef HAVE_MBA_3D
      <<"yes\n"
#else
      <<"no\n"
#endif
      <<"CGAL support\t\t\t"
#ifdef HAVE_CGAL
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
      <<"Hypre support\t\t\t"
#ifdef HAVE_HYPRE
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
      <<"Python support\t\t\t"
#ifdef HAVE_PYTHON
      <<"yes\n"
#else
      <<"no\n"
#endif
      <<"Numpy support\t\t\t"
#ifdef HAVE_NUMPY
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
      <<"Zoltan support\t\t\t"
#ifdef HAVE_ZOLTAN
      <<"yes\n"
#else
      <<"no\n"
#endif
      ;
  stream.flush();
  return;
}

void print_environment(){

  const char *relevant_variables[]={"FLUIDITY_INTEGER_SIZE", "FLUIDITY_REAL_SIZE",
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

void usage(char *cmd){
  print_version();
  cerr<<"\n\nUsage: "<<cmd<<" [options ...] [simulation-file]\n"
      <<"\nOptions (NOTE: Long options are not available on AIX.):\n"
      <<" -a <mode>, --adjoint <mode>\n\tRun the full adjoint model. Possible modes are {forward, initial, bc, sensitivity, parameter,getobservation}\n"
      <<" -h, --help\n\tHelp! Prints this message.\n"
      <<" -k, --vtk\n\tDumps out VTK files (.vtu) rather than the usual .d files.\n" 
      <<" -l, --log\n\tCreate log file for each process (useful for non-interactive testing)."
      <<" Sets default value for -v to 2.\n"
      <<" -r <mode>, --reduced <mode>\n\tRun the POD reduced model. Possible modes are {forward, initial, bc, sensitivity, parameter,getobservation}\n"
      <<" -v <level>, --verbose\n\tVerbose output to stdout, default level 0\n"
      <<" -V, --version\n\tVersion\n"
      <<" -d, --pseudo2d\n\tSpecify squashed dimension for pseudo2d domains\n"
      <<" -p, --petsc\n\tPetsc Options (e.g. -p -info,-geopressure_ksp_type,cg)\n"
      <<" -x, --xml\n\tRun from a new xml options file\n";
  return;
}

void ParseArguments(int argc, char** argv){

#ifndef _AIX
  struct option longOptions[] = {
    {"adjoint", 1, 0, 'a'},
    {"help", 0, 0, 'h'},
    {"vtk", 0, 0, 'k'},
    {"log", 0, 0, 'l'},
    {"verbose", optional_argument, 0, 'v'},
    {"version", 0, 0, 'V'},
    {"pseudo2d", 0, 0, 'd'},
    {"petscopts", 0, 0, 'p'},
    {"reduced", 1, 0, 'r'},
    {"xml", 1, 0, 'x'},
    {0, 0, 0, 0}
  };
#endif
  int optionIndex = 0;
  int verbosity = 0;
  int c;

  // set opterr to nonzero to make getopt print error messages 
  opterr=1;

  while (true){
#ifndef _AIX
    c = getopt_long(argc, argv, "a:hklr:v::Sd:p:Vx:", longOptions, &optionIndex);
#else
    c = getopt(argc, argv, "a:hklr:v::Sd:p:Vx:");
#endif
    if (c == -1) break;

    OptionError stat;
    switch (c){
    case 'a':
      // Add this to the options tree
      stat = set_option(string("/model/fluids/adjoint/") + string(optarg), "enable");
      assert(stat == SPUD_NO_ERROR or stat == SPUD_NEW_KEY_WARNING);
      break;

    case 'h':
      fl_command_line_options["help"] = "";
      break;

    case 'k':
      stat = set_option("io/dump_format", "vtk");
      assert(stat == SPUD_NO_ERROR or stat == SPUD_NEW_KEY_WARNING);
      break;
      
    case 'l':
      fl_command_line_options["log"] = "";
      break;        
      
    case 'r':
#if defined(HAVE_LIBARPACK)
      // Add this to the options tree
      stat = set_option(string("/model/fluids/reduced/") + string(optarg), "enable");
      assert(stat == SPUD_NO_ERROR or stat == SPUD_NEW_KEY_WARNING);   
#else
      cerr<<"ERROR: missing ARPACK support\n";
      usage(argv[0]);
      exit(-1);
#endif
      break;

    case 'v':
      fl_command_line_options["verbose"] = (optarg == NULL) ? "1" : optarg;
      break;  
      
    case 'V':
      fl_command_line_options["version"] = "";
      break;
    
    case 'd':
      fl_command_line_options["pseudo2d"] = optarg;
      break;
    
    case 'p':
      fl_command_line_options["petscopts"] = optarg;
      break;

    case 'x':
      fl_command_line_options["xml"] = optarg;
      break;

    case '?':
      // missing argument only returns ':' if the option string starts with ':'
      // but this seems to stop the printing of error messages by getopt?
      cerr<<"ERROR: unknown option or missing argument\n";
      usage(argv[0]);
      exit(-1);
    case ':':
      cerr<<"ERROR: missing argument\n";
      usage(argv[0]);
      exit(-1);
    default:
      // unexpected:
      cerr<<"ERROR: getopt returned unrecognized character code\n";
      exit(-1);
    }
  }
  
  // Help?
  if(fl_command_line_options.count("help")){
    usage(argv[0]);
    exit(-1);
  }
  
  // Version?
  if(fl_command_line_options.count("version")){
    print_version();
    exit(-1);
  }
  
  // Verbose?
  {    
    int MyRank = 0;
#ifdef HAVE_MPI
    if(MPI::Is_initialized()){
      MyRank = MPI::COMM_WORLD.Get_rank();
    }
#endif
 
    if(fl_command_line_options.count("verbose") == 0){
      verbosity = 0;
    }else{
      verbosity = atoi(fl_command_line_options["verbose"].c_str());
    }
    set_global_debug_level_fc(&verbosity);
  }
  
  // Pseudo2d?
  if(fl_command_line_options.count("pseudo2d")){
  int val;
  val = atoi(fl_command_line_options["pseudo2d"].c_str());
  set_pseudo2d_domain_fc(&val);
  }

  // What to do with stdout/stderr?
  if(fl_command_line_options.count("log")){
    ostringstream debug_file, err_file;
    debug_file << "fluidity.log";
    err_file << "fluidity.err";
#ifdef HAVE_MPI
    if(MPI::Is_initialized()){
      int MyRank = MPI::COMM_WORLD.Get_rank();
      debug_file << "-" << MyRank;
      err_file << "-" << MyRank;
    }
#endif
  
    if(freopen(debug_file.str().c_str(), "w", stdout) == NULL)
      perror("failed to redirect stdio for debugging"); 
        
    if(freopen(err_file.str().c_str(), "w", stderr) == NULL)
      perror("failed to redirect stderr for debugging");
  }
  
  // Find the filename if one is specified. Assuming that final
  // argument is simulation filename if it does not correspond to a
  // known option.
  if (fl_command_line_options.count("xml") == 0)
  {
    if(argc > optind + 1)
    {
      fl_command_line_options["xml"] = argv[optind + 1];
    }
    else if(argc == optind + 1)
    {
      fl_command_line_options["xml"] = argv[optind];
    }
    else if(fl_command_line_options.count("xml") == 0)
    {
      cerr << "ERROR: Unrecognized arguments!" << endl;
      usage(argv[0]);
      exit(-1);
    }
  }
  
  load_options(fl_command_line_options["xml"]);
  if(!have_option("/simulation_name")){
    cerr<<"ERROR: failed to find simulation name after loading options file\n";
    cerr<<"Note that running from gem files is no longer supported.\n";
    exit(-1);
  }
  OptionError stat = get_option("/simulation_name", fl_command_line_options["simulation_name"]);
  assert(stat == SPUD_NO_ERROR);
  
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

  // Environmental stuff -- this needs to me moved out of here and
  // into populate state.
  if(have_option("/timestepping/current_time/time_units")){
    string option;
    get_option("/timestepping/current_time/time_units/date", option);

    FluxesReader_global.SetSimulationTimeUnits(option.c_str());
    ClimateReader_global.SetSimulationTimeUnits(option.c_str());
    NEMOReader_v2_global.SetSimulationTimeUnits(option.c_str());
  }
  if(have_option("/environmental_data/climatology/file_name")){
    string option;
    get_option("/environmental_data/climatology/file_name", option);

    ClimateReader_global.SetClimatology(option);
  }
  if(have_option("/ocean_forcing/input_file")) {
    string option;
    get_option("/ocean_forcing/input_file/file_name", option);

    FluxesReader_global.RegisterDataFile(option);
#ifdef DDEBUG
    //FluxesReader_global.VerboseOn();
#endif
    // field from NetCDF file          Index |   Physical meaning
    FluxesReader_global.AddFieldOfInterest("10u");  //  0   | 10 metre U wind component
    FluxesReader_global.AddFieldOfInterest("10v");  //  1   | 10 metre V wind component
    FluxesReader_global.AddFieldOfInterest("ssrd"); //  2   | Surface solar radiation
    FluxesReader_global.AddFieldOfInterest("strd"); //  3   | Surface thermal radiation 
    FluxesReader_global.AddFieldOfInterest("ro");   //  4   | Runoff
    FluxesReader_global.AddFieldOfInterest("tp");   //  5   | Total precipitation
    FluxesReader_global.AddFieldOfInterest("2d");   //  6   | Dew point temp at 2m
    FluxesReader_global.AddFieldOfInterest("2t");   //  7   | Air temp at 2m 
    FluxesReader_global.AddFieldOfInterest("msl");  //  8   | Mean sea level pressure 
  }

  if(have_option("/ocean_forcing/external_data_boundary_conditions")) {
    string option;
    get_option("/ocean_forcing/external_data_boundary_conditions/input_file/file_name", option);

    if(verbosity >= 3)
      cout << "Registering external data forcing file: " << option << endl;

    NEMOReader_v2_global.RegisterDataFile(option);

    NEMOReader_v2_global.AddFieldOfInterest("temperature");  //  0   | Sea temperature
    NEMOReader_v2_global.AddFieldOfInterest("salinity");     //  1   | Salinity
    NEMOReader_v2_global.AddFieldOfInterest("u");            //  2   | Azimuthal velocity
    NEMOReader_v2_global.AddFieldOfInterest("v");            //  3   | Meridional velocity
    NEMOReader_v2_global.AddFieldOfInterest("ssh");          //  4   | Sea surface height
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
    for(int i = 1; i < petscargc; i++)
      petscargv[i] = new char[petscopts[i-1].size()+1];
    strncpy(petscargv[0], argv[0], strlen(argv[0]) + 1);
    for(int i = 1; i < petscargc; i++)
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
  
  for(int i = 0; i < petscargc; i++)
    delete [] petscargv[i];
  delete [] petscargv;
#endif
}

