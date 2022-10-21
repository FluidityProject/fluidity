/* h5pToGNUplot.cc
   Antino Kim
   This utility will output a GNU Plot input file accodring to the flags provided from the command line.
   The parser was imported from the example of h5dump utility with slight modifications.
*/

#include <stdio.h>
#include <stdlib.h>
#include <cctype>
#include <string.h>
#include <hdf5.h>
#include "H5Part.h"

#define MAX_LEN 100

/* Function headers */
int get_option(int argc, const char **argv, const char *opts, const struct long_options *l_opts);
static void print_help();
static void variable_assign(int argc, const char *argv[]);

/* Global variables */
static char* var_1      = NULL;
static char* var_2      = NULL;
static char* input_name      = NULL;
static char* output_name      = NULL;
static char* timestep      = NULL;
static char* start      = NULL;
static char* npoints      = NULL;
static int print_all = 0;

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
static const char *s_opts = "h1:2:i:o:n:t:s:";

/* List of options in full words */
static struct long_options l_opts[] =
{
    { "help", no_arg, 'h' },         // Print help page
    { "1var", require_arg, '1' },         // Takes first variable parameter
    { "2var", require_arg, '2' },        // Takes second variable parameter
    { "input", require_arg, 'i' },        // Takes input file name
    { "output", require_arg, 'o' },   // Takes output file name (without this flag, the program will print to stdout)
    { "number", require_arg, 'n' }, // Sets number of output points
    { "start", require_arg, 's' },  // Sets the starting particle index
    { "timestep", require_arg, 't' },  // Sets the timestep
    { NULL, 0, '\0' }
};


/************************************************************************************
***********************************  FUNCTIONS  *************************************
*************************************************************************************/


/* get_option is the parsing function that was majorly ported from h5dump utility */
int get_option(int argc, const char **argv, const char *opts, const struct long_options *l_opts)
{
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
        register char *cp;    /* pointer into current token */

        /* short command line option */
        opt_opt = argv[opt_ind][sp];

        if (opt_opt == ':' || (cp = strchr(opts, opt_opt)) == 0)
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
          case '1': // Print file attributes
            var_1 = strdup(opt_arg);
            break;
          case '2': // Print step attributes & values for time step n
            var_2 = strdup(opt_arg);
            break;
          case 'i': // Print data sets names & values for time step n
            input_name = strdup(opt_arg);
            break;
          case 'o': // Print number of steps
            output_name = strdup(opt_arg);
            break;
          case 'n': // Print shorter version without the values
            npoints = strdup(opt_arg);
            break;
          case 's': // Print shorter version without the values
            start = strdup(opt_arg);
            break;
          case 't': // Print shorter version without the values
            timestep = strdup(opt_arg);
            if(atoi(timestep)==-1)
            {
               print_all = 1;
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
   fprintf(stdout, "\nusage: h5pToGNUplot -t TIMESTEP -1 VARIABLE#1 -2 VARIABLE#2 -i INPUTFILE [OPTIONAL_FLAGS]\n");
   fprintf(stdout, "\n");
   fprintf(stdout, "  FLAGS\n");
   fprintf(stdout, "   -h, --help                 Print help page\n");
   fprintf(stdout, "   -1 par, --1var par         (REQUIRED) Takes first variable parameter to \"par\"\n");
   fprintf(stdout, "   -2 par, --2var par         (REQUIRED) Takes second variable parameter to \"par\"\n");
   fprintf(stdout, "   -i file, --input file      (REQUIRED) Takes input file name to \"file\"\n");
   fprintf(stdout, "   -t step, --timestep step   (REQUIRED) Sets the timestep to \"step\" (Value -1 will result in dumping values of all timesteps.)\n");
   fprintf(stdout, "   -o file, --output file     (OPTIONAL) Takes output file name to \"file\" (without this flag, the program will print to stdout)\n");
   fprintf(stdout, "   -n num, --number num       (OPTIONAL) Sets number of output points to \"num\"\n");
   fprintf(stdout, "   -s idx, --start idx        (OPTIONAL) Sets the starting particle index to \"idx\"\n");
   fprintf(stdout, "\n");
   fprintf(stdout, "  Examples:\n");
   fprintf(stdout, "\n");
   fprintf(stdout, "  1) Create GNU plot file output.txt from sample.h5part by ploting x vs px for timestep 54\n");
   fprintf(stdout, "\n");
   fprintf(stdout, "        h5pToGNUplot -t 54 -1 x -2 px -i sample.h5part -o output.txt\n");
   fprintf(stdout, "\n");
   fprintf(stdout, "  2) Create GNU plot file output.txt from sample.h5part by ploting x vs px for timestep 54 \n     using 1200 points from particle index 76\n");
   fprintf(stdout, "\n");
   fprintf(stdout, "        h5pToGNUplot -t 54 -1 x -2 px -i sample.h5part -o output.txt -s 76 -n 1200\n");
   fprintf(stdout, "\n");
}

int main(int argc, const char *argv[])
{
   /* Numerous variables */
   H5PartFile           *h5file = NULL;
   char data_name[MAX_LEN];
   int i, j, start_indx = 0;
   h5part_int64_t type;
   hid_t type_1;
   hid_t type_2;
   int num_dataset;
   h5part_int64_t nparticles;
   int stop = 0;
   void* value_1 = NULL;
   void* value_2 = NULL;
   FILE * outFile;
   int ntime_step = 0;

   variable_assign(argc, argv);

   if(input_name == NULL) 
   {
       fprintf(stdout, "missing file name\n");
       print_help();
       exit(1);
   }

  if(var_1 == NULL || var_2 == NULL || timestep == NULL)
   {
       fprintf(stdout, "missing required parameters\n");
       print_help();
       exit(1);
   }

   h5file = H5PartOpenFile(input_name, H5PART_READ);

   if( h5file == NULL )
   {
        fprintf(stdout, "unable to open file %s\n", input_name);
        print_help();
        exit(1);
   }

   j = atoi(timestep);
   ntime_step = j+1;

   if(print_all)
   {
      j=0;
      ntime_step = H5PartGetNumSteps(h5file);
   }

   for (; j<ntime_step; j++)
   {
      H5PartSetStep(h5file,j);
      num_dataset = H5PartGetNumDatasets(h5file);

      for(i=0; i<num_dataset; i++)
      {
         H5PartGetDatasetInfo(h5file, i, data_name, MAX_LEN, &type, &nparticles);
         if(strcmp(data_name, var_1))
         {
            type_1=type;
         }
         if(strcmp(data_name, var_2))
         {
            type_2=type;
         }
      }

      if (type_1 == H5T_NATIVE_INT64)
      {
         value_1 = (h5part_int64_t*)malloc(sizeof(h5part_int64_t)*nparticles);
         H5PartReadDataInt64(h5file, var_1, (h5part_int64_t*)value_1);
      }
      else if (type_1 == H5T_NATIVE_DOUBLE)
      {
         value_1 = (double*)malloc(sizeof(double)*nparticles);
         H5PartReadDataFloat64(h5file, var_1, (double*)value_1);
      }
      else
      {
         fprintf(stdout, "Dataset Type is UNKNOWN for %s. (Check the variable name.)\n", var_1);
         exit(1);
      }

      if (type_2 == H5T_NATIVE_INT64)
      {
         value_2 = (h5part_int64_t*)malloc(sizeof(h5part_int64_t)*nparticles);
         H5PartReadDataInt64(h5file, var_2, (h5part_int64_t*)value_2);
      }
      else if (type_2 == H5T_NATIVE_DOUBLE)
      {
         value_2 = (double*)malloc(sizeof(double)*nparticles);
         H5PartReadDataFloat64(h5file, var_2, (double*)value_2);
      }
      else
      {
         fprintf(stdout, "Dataset Type is UNKNOWN for %s. (Check the variable name.)\n", var_2);
         exit(1);
      }

     if(start != NULL)
     {
       start_indx = atoi(start);
     }

     if(npoints != NULL)
     {
       stop = atoi(npoints) + start_indx;
     }
     else
     {
       stop = nparticles;
     }

     if(stop > nparticles)
     {
       stop = nparticles;
     }


     if(output_name != NULL)
     {
       outFile = fopen (output_name,"w");

       for(i=start_indx; i<stop; i++)
       {
          if (type_1 == H5T_NATIVE_INT64)
            fprintf(outFile, "%lld", ((long long*)value_1)[i]);
          else if (type_1 == H5T_NATIVE_DOUBLE)
            fprintf(outFile, "%lf", ((double*)value_1)[i]);

          if (type_2 == H5T_NATIVE_INT64)
            fprintf(outFile, "\t%lld", ((long long*)value_2)[i]);
          else if (type_2 == H5T_NATIVE_DOUBLE)
            fprintf(outFile, "\t%lf", ((double*)value_2)[i]);

         fprintf(outFile, "\n");
       }
       fprintf(outFile, "\n");
       fclose (outFile);
     }
     else
     {
       for(i=start_indx; i<stop; i++)
       {
          if (type_1 == H5T_NATIVE_INT64)
            fprintf(stdout, "%lld", ((long long*)value_1)[i]);
          else if (type_1 == H5T_NATIVE_DOUBLE)
            fprintf(stdout, "%lf", ((double*)value_1)[i]);

          if (type_2 == H5T_NATIVE_INT64)
            fprintf(stdout, "\t%lld", ((long long*)value_2)[i]);
          else if (type_2 == H5T_NATIVE_DOUBLE)
            fprintf(stdout, "\t%lf", ((double*)value_2)[i]);

         fprintf(stdout, "\n");
       }
       fprintf(stdout, "\n");
     }

     free(value_1);
     free(value_2);
  }

  H5PartCloseFile(h5file);
  return 0;
}
