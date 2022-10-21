/* h5PartSurfaceToVtk.cc
   Andreas Adelmann & Chuan Wang
   2010-2011
*/

#include <hdf5.h>
#include "H5hut.h"

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <cctype>
#include <string.h>

#include <set>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cassert>

using namespace std;

#define MAX_LEN 100

/* Function headers */
int get_option(int argc, const char **argv, const char *opts, const struct long_options *l_opts);

static void print_help();

static void variable_assign(int argc, const char *argv[]);

/* Global variables */

static char* input_name      = NULL;
static char* mytemp          = NULL;
static char* geo_name        = NULL;
static char* output_name     = NULL;

static bool flg_alive = false;
static double z_pos = 0.0; 
static int print_all = 0;
static int step_ini = 1;
static int dump_freq = 1;


/* `get_option' variables */
int         opt_err = 1;    /*get_option prints errors if this is on */
int         opt_ind = 1;    /*token pointer                          */
const char *opt_arg = NULL;        /*flag argument (or value)               */

/* indication whether the flag (option) requires an argument or not */
enum {
  no_arg = 0,         /* doesn't take an argument     */
  require_arg,        /* requires an argument	        */
};

/* struct for flags (options) */
typedef struct long_options
{
  const char  *name;          /* name of the long option              */
  int          has_arg;       /* whether we should look for an arg    */
  char         shortval;      /* the shortname equivalent of long arg
			       * this gets returned from get_option   */
} long_options;

/* List of options in single characters */
static const char *s_opts = "h1:2:i:g:o:f:d:";

/* List of options in full words */
static struct long_options l_opts[] =
  {
    { "help", no_arg, 'h' },         // Print help page
    { "input", require_arg, 'i' },        // Takes input file name
    { "geometry", require_arg, 'g' },        // Takes input geometry file name
    { "output", require_arg, 'o' },   // Takes output file name (without this flag, the program will print to stdout)
    { "first", require_arg, 'f' },   // first step
    { "dumpfreq", require_arg, 'd' },   // dump frequency of the
    { NULL, 0, '\0' }
  };

/************************************************************************************
 ***********************************  FUNCTIONS  *************************************
 *************************************************************************************/


string convert2Int(int number) {
  stringstream ss;
  ss <<   setw(5) << setfill('0') <<  number; 
  return ss.str(); 
}


/* get_option is the parsing function that was majorly ported from h5dump utility */
int get_option(int argc, const char **argv, const char *opts, const struct long_options *l_opts) {
  static int sp = 1;    /* character index in current token */
  int opt_opt = '?';    /* option character passed back to user */

  if (sp == 1) 
    {
      /* check for more flag-like tokens */
      if (opt_ind >= argc || argv[opt_ind][0] != '-' || argv[opt_ind][1] == '\0') 
	{
	  return EOF;
	}
      else if (strcmp(argv[opt_ind], "--") == 0)
	{
	  opt_ind++;
	  return EOF;
	}
    }

  if (sp == 1 && argv[opt_ind][0] == '-' && argv[opt_ind][1] == '-') 
    {
      /* long command line option */
      const char *arg = &argv[opt_ind][2];
      int i;

      for (i = 0; l_opts && l_opts[i].name; i++)
	{
	  size_t len = strlen(l_opts[i].name);

	  if (strncmp(arg, l_opts[i].name, len) == 0)
	    {
	      /* we've found a matching long command line flag */
	      opt_opt = l_opts[i].shortval;

	      if (l_opts[i].has_arg != no_arg)
		{
		  if (arg[len] == '=')
		    {
		      opt_arg = &arg[len + 1];
		    }
		  else if (opt_ind < (argc - 1) && argv[opt_ind + 1][0] != '-')
		    {
		      opt_arg = argv[++opt_ind];
		    }
		  else if (l_opts[i].has_arg == require_arg)
		    {
		      if (opt_err)
			fprintf(stderr, "%s: option required for \"--%s\" flag\n", argv[0], arg);

		      opt_opt = '?';
		    }
		}
	      else
		{
		  if (arg[len] == '=')
		    {
		      if (opt_err)
			fprintf(stderr, "%s: no option required for \"%s\" flag\n", argv[0], arg);

		      opt_opt = '?';
		    }

		  opt_arg = NULL;
		}

	      break;
	    }
	}

      if (l_opts[i].name == NULL)
	{
	  /* exhausted all of the l_opts we have and still didn't match */
	  if (opt_err)
	    fprintf(stderr, "%s: unknown option \"%s\"\n", argv[0], arg);

	  opt_opt = '?';
	}

      opt_ind++;
      sp = 1;
    }
  else
    {
      register const char *cp;    /* pointer into current token */

      /* short command line option */
      opt_opt = argv[opt_ind][sp];

      if (opt_opt == ':' || (cp = (char *)strchr(opts, opt_opt)) == 0)
	{

	  if (opt_err)
	    fprintf(stderr, "%s: unknown option \"%c\"\n", argv[0], opt_opt);
	  /* if no chars left in this token, move to next token */
	  if (argv[opt_ind][++sp] == '\0')
	    {
	      opt_ind++;
	      sp = 1;
	    }

	  return '?';
	}

      if (*++cp == ':')
	{

	  /* if a value is expected, get it */
	  if (argv[opt_ind][sp + 1] != '\0')
	    {
	      /* flag value is rest of current token */
	      opt_arg = &argv[opt_ind++][sp + 1];
	    }
	  else if (++opt_ind >= argc)
	    {
	      if (opt_err)
		{
		  fprintf(stderr, "%s: value expected for option \"%c\"\n", argv[0], opt_opt);
		}
	      opt_opt = '?';
	    }
	  else
	    {
	      /* flag value is next token */
	      opt_arg = argv[opt_ind++];
	    }

	  sp = 1;
	}
      else 
	{
	  /* set up to look at next char in token, next time */
	  if (argv[opt_ind][++sp] == '\0')
	    {
	      /* no more in current token, so setup next token */
	      opt_ind++;
	      sp = 1;
	    }

	  opt_arg = NULL;
	}
    }

  /* return the current flag character found */
  return opt_opt;
}

/* Assigns functions according to the parsed result */
static void variable_assign(int argc, const char *argv[])
{
  int option;

  /* set options according to the command line */
  while ((option = get_option(argc, argv, s_opts, l_opts)) != EOF)
    {
      switch ((char)option)
	{
	case 'h': // Print help page
	  print_help();
	  exit(1);
	case 'o': // Print number of steps
	  output_name = strdup(opt_arg);
	  break;
	case 'i': // Print shorter version without the values
	  input_name = strdup(opt_arg);
	  break;
	case 'g': // 
	  { geo_name = strdup(opt_arg);
                        
	  }
	  break;  
	case 'f': // 
	  { 	
	    mytemp = strdup(opt_arg);
	    step_ini = (int)strtod(mytemp,NULL);
	    fprintf(stdout,"string=%s, %s \n",opt_arg, mytemp);
	    fprintf(stdout,"step_ini=%i \n",step_ini);                        
	  }
	  break;  
	case 'd': // 
	  { 	
	    fprintf(stdout,"string=%s \n",opt_arg);
	    mytemp = strdup(opt_arg);
	    dump_freq = (int)strtod(mytemp,NULL);
	    fprintf(stdout,"string=%s \n",opt_arg);
	    fprintf(stdout,"dump_freq=%i \n",dump_freq);
	  }
	  break;  
	default:
	  print_help();
	  exit(1);
	}
    }
}

/* For printing help page */
static void print_help()
{
  fflush(stdout);
  fprintf(stdout, "\nusage: h5PartSurfaceToVtk  -i INPUTFILE -g GEOMETRYFILE -o OUTPUTFILE\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "  FLAGS\n");
  fprintf(stdout, "   -h, --help                 Print help page\n");
  fprintf(stdout, "   -i file, --input file      (REQUIRED) Takes input base file name to \"file\" (extension h5 is assumed \n");
  fprintf(stdout, "   -g file, --input geometry file      (REQUIRED) Takes input base file name to \"file\" (extension h5 is assumed \n");
  fprintf(stdout, "   -o file, --output file     (REQUIRED) Takes output base file name to \"file\" (extension vtk is added)\n");
    
  fprintf(stdout, "\n");
  fprintf(stdout, "  Examples first dump step 100 and dump frequency 100:\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "        /h5SurfaceVtk -i Cavity_Mag_Surface -g ../10 -o p- -f 100 -d 100 \n");
  fprintf(stdout, "\n");
    
}


int main(int argc, const char *argv[])
{
  h5_file_t *h5file = NULL;
    
  std::ofstream of, ofalive, ofenergy, ofpnum;
   
  int j;

  int num_dataset;

  int ntime_step = 0;

  variable_assign(argc, argv);

  if(input_name == NULL)  {
    fprintf(stdout, "missing input file name\n");
    print_help();
    exit(1);
  }
   
  if(output_name == NULL) {
    fprintf(stdout, "missing output file name\n");
    print_help();
    exit(1);
  }
   
  string ifn = string(input_name) + string(".h5");
   
  h5file = H5OpenFile(ifn.c_str(), H5_O_WRONLY, 0);
     
  if( h5file == NULL ) {
    fprintf(stdout, "unable to open file %s\n", input_name);
    print_help();
    exit(1);
  }

  string gfn = string(geo_name) + string(".h5");
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS); //Property list identifier
  //    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  // Create a new file collectively and release property list identifier.
  hid_t file_id = H5Fopen(gfn.c_str(), H5F_ACC_RDONLY, plist_id);
  assert(file_id >= 0);
  H5Pclose(plist_id);
  int     numbfaces_global_m;
  /////////////////////////////////////////////
    //   Read dataset "surface" from .h5 file
    ////////////////////////////////////////////

    hsize_t dimsf[2];//dataset dimensions
    herr_t  status;

    hid_t dset_id = H5Dopen(file_id, "/surface",H5P_DEFAULT);
    assert(dset_id >= 0);
    // Create the dataspace for the dataset.
    hid_t x = H5Dget_space(dset_id);
    H5Sget_simple_extent_dims(x, dimsf, NULL);
    hid_t filespace = H5Screate_simple(2, dimsf, NULL);

    numbfaces_global_m = dimsf[0];
    int *allbfaces_m = (int*)malloc(numbfaces_global_m * 4*sizeof(int));

    int *bfaces_idx_m = (int*)malloc(numbfaces_global_m*sizeof(int));

    //Tribarycent_m = new Vector_t[numbfaces_global_m];
    // Create property list for collective dataset write.
    hid_t plist_id2 = H5Pcreate(H5P_DATASET_XFER);
    //H5Pset_dxpl_mpio(plist_id2, H5FD_MPIO_COLLECTIVE);
    status = H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, filespace, plist_id2, allbfaces_m);
    assert(status >= 0);
    // Store local boundary faces, discard others
    int nof_sym_m = 0; // Count number of symmetry planes.
    int const nogeosym_flag = 0; // heronion outputs a index, case 0 indicates no symmetry plane, 1 for (x,y) sym, 2 for (y,z), 3 for (z,x),none 0 means the current triangle is not on a real surface of geometry but on a symmetry plane.

   
    for(int i = 0; i < numbfaces_global_m; i ++) {
      bfaces_idx_m[i] = i;
      if(allbfaces_m[4 * i] > nogeosym_flag) {
	nof_sym_m += 1;
	if(i < numbfaces_global_m - 1) {
	  for(int j = 0; j < numbfaces_global_m - i; j ++) {
	    allbfaces_m[4*(i+j)] = allbfaces_m[4*(i+j+1)];
	    allbfaces_m[4*(i+j)+1] = allbfaces_m[4*(i+j+1)+1];
	    allbfaces_m[4*(i+j)+2] = allbfaces_m[4*(i+j+1)+2];
	    allbfaces_m[4*(i+j)+3] = allbfaces_m[4*(i+j+1)+3];
	  }
	} else
	  numbfaces_global_m = numbfaces_global_m - 1;
      }
    }
    H5Dclose(dset_id);
    H5Sclose(filespace);
    ////////////////////////////////
    //  Also read dataset "coords"
    ////////////////////////////////
    hsize_t dimsf_c[2];
    herr_t status_c;
    hid_t dset_id_c = H5Dopen(file_id, "/coords",H5P_DEFAULT);
    assert(dset_id_c >= 0);

    // Create the dataspace for the dataset.
    hid_t cp_space = H5Dget_space(dset_id_c);
    H5Sget_simple_extent_dims(cp_space, dimsf_c, NULL);
    hid_t filespace_c = H5Screate_simple(2, dimsf_c, NULL);

    int numpoints_global_m = dimsf_c[0];
    double *point_coords = (double*)malloc(3 * numpoints_global_m*sizeof(double));
    hid_t plist_id3 = H5Pcreate(H5P_DATASET_XFER);
    //H5Pset_dxpl_mpio(plist_id3, H5FD_MPIO_COLLECTIVE);
    status_c = H5Dread(dset_id_c, H5T_NATIVE_DOUBLE, H5S_ALL, filespace_c, plist_id3, point_coords);
    assert(status_c >= 0);
    

    
    
    double **geo3Dcoords_m;//=  malloc(sizeof(double)*numpoints_global_m*3);
    geo3Dcoords_m = (double**)malloc(sizeof(double)*numpoints_global_m);
    for (int i=0;i<numpoints_global_m;i++) {
      geo3Dcoords_m[i] = (double*)malloc(sizeof(double)*3);
    }
    
    for(int i = 0; i < numpoints_global_m; i++) {
      geo3Dcoords_m[i][0] = point_coords[3*i];
      geo3Dcoords_m[i][1] = point_coords[3*i+1];
      geo3Dcoords_m[i][2] = point_coords[3*i+2];
    }
   
    
    
    // Close HDF5 stuff
    H5Dclose(dset_id_c);
    H5Sclose(filespace_c);


    ntime_step = H5GetNumSteps(h5file);
    cout<<"step number: "<<ntime_step<<endl;
    double sey_num;//sey total number
    int remained_num;//
    for (size_t j=1; j<=ntime_step; j++) {
      set<h5_int64_t> idSet;
      size_t step_local = step_ini + (j-1)*dump_freq;
      H5SetStep(h5file,step_local);
      num_dataset = H5PartGetNumDatasets(h5file);
      h5_int64_t triNum = H5PartGetNumParticles(h5file);
      double    qi = 0.0;
      H5ReadStepAttribFloat64(h5file,"qi",&qi);
      cout << "Working on timestep " << step_local  << endl;
		
      h5_float64_t* PrimaryLoss = (h5_float64_t*)malloc(sizeof(h5_float64_t)*triNum);
      H5PartReadDataFloat64(h5file, "PrimaryLoss", PrimaryLoss);
	
      h5_float64_t* SecondaryLoss = (h5_float64_t*)malloc(sizeof(h5_float64_t)*triNum);
      H5PartReadDataFloat64(h5file, "SecondaryLoss", SecondaryLoss);
	
      h5_float64_t* FNEmissionLoss = (h5_float64_t*)malloc(sizeof(h5_float64_t)*triNum);
      H5PartReadDataFloat64(h5file, "FNEmissionLoss", FNEmissionLoss);
	
      h5_int64_t* TriID = (h5_int64_t*)malloc(sizeof(h5_int64_t)*triNum);
      H5PartReadDataInt64(h5file, "TriangleID", TriID);
	
      string ffn = string("vtk/") + string(output_name) + string("Surface-") + convert2Int(j) + string(".vtk");
	
      of.open(ffn.c_str());
      assert(of.is_open());
      of.precision(6);
      of << "# vtk DataFile Version 2.0" << endl;  // first line of a vtk file
      of << "generated using DataSink::writeGeoToVtk" << endl;   // header
      of << "ASCII" << endl << endl;   // file format
      of << "DATASET UNSTRUCTURED_GRID" << endl;  // dataset structrue
      of << "POINTS " << numpoints_global_m << " float" << endl;   // data num and type
      for(int i = 0; i < numpoints_global_m ; i ++)
	of << geo3Dcoords_m[i][0] << " " << geo3Dcoords_m[i][1] << " " << geo3Dcoords_m[i][2] << endl;
      of << endl;

      of << "CELLS " << numbfaces_global_m << " " << 4 * numbfaces_global_m << endl;
      for(int i = 0; i < numbfaces_global_m ; i ++)
	of << "3 " << allbfaces_m[4*i+1] << " " << allbfaces_m[4*i+2] << " " << allbfaces_m[4*i+3] << endl;
      of << "CELL_TYPES " << numbfaces_global_m << endl;
      for(int i = 0; i < numbfaces_global_m ; i ++)
	of << "5" << endl;
      of << "CELL_DATA " << numbfaces_global_m << endl;
      of << "SCALARS " << "PrimaryLoss" << " float " << "1" << endl;
      of << "LOOKUP_TABLE " << "default" << endl ;
      for(int i = 0; i < numbfaces_global_m ; i ++)
	of << fabs(PrimaryLoss[i]/qi) << endl;
      of << "SCALARS " << "SecondaryLoss" << " float " << "1" << endl;
      of << "LOOKUP_TABLE " << "default" << endl ;
      for(int i = 0; i < numbfaces_global_m ; i ++)
	of << fabs(SecondaryLoss[i]/qi) << endl;
      of << "SCALARS " << "FNEmissionLoss" << " float " << "1" << endl;
      of << "LOOKUP_TABLE " << "default" << endl ;
      for(int i = 0; i < numbfaces_global_m ; i ++)
	of << fabs(FNEmissionLoss[i]/qi) << endl;
      of << "SCALARS " << "TotalLoss" << " float " << "1" << endl;
      of << "LOOKUP_TABLE " << "default" << endl ;
      for(int i = 0; i < numbfaces_global_m ; i ++)
	of << fabs((FNEmissionLoss[i]+SecondaryLoss[i]+PrimaryLoss[i])/qi) << endl;
      of << endl;
	
      of.close();
	   
      //cout <<"ptype size = "<< sizeof(ptype) << endl;
	
      free(PrimaryLoss);
      free(SecondaryLoss);
      free(FNEmissionLoss);
      free(TriID);
    }
    free(point_coords);
    for (int i=0;i<numpoints_global_m;i++) {
      free(geo3Dcoords_m[i]);
    }
    free(geo3Dcoords_m);
    free(allbfaces_m);
    free(bfaces_idx_m);

    H5CloseFile(h5file);   
    return 0;
}
