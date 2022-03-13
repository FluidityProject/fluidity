/* h5pAttrib.cc
   Antino Kim
   This utility will output information on h5part files accodring to the flags provided from the command line.
   The parser was imported from the example of h5dump utility with slight modifications.
*/

#include <stdio.h>
#include <stdlib.h>
#include <cctype>
#include <string.h>

#include "H5hut.h"

#define MAX_LEN 100

/* Function headers */
int get_option(int argc, const char **argv, const char *opts, const struct long_options *l_opts);
static void print_help();
static void free_handler(struct arg_handler *hand, int len);
static void print_all(h5_file_t* file);
static void print_nstep(h5_file_t* file, char * garbage);
static void print_file_attributes(h5_file_t* file, char * garbage);
static void print_step_attributes(h5_file_t* file, char *attr);
static void print_dataset(h5_file_t* file, char *attr);
static struct arg_handler* function_assign(int argc, const char *argv[]);

/* Global variables */
static int         display_all       = true;
static int         print_header      = false;
static const char* global_fname      = NULL;

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
static const char *s_opts = "hnAHa:d:";

/* List of options in full words */
static struct long_options l_opts[] =
{
    { "help", no_arg, 'h' },         // Print help page
    { "nstep", no_arg, 'n' },        // Print number of steps
    { "fileA", no_arg, 'A' },        // Print file attributes
    { "stepA", require_arg, 'a' },   // Print step attributes & values for time step n
    { "dataset", require_arg, 'd' }, // Print data sets names & values for time step n
    { "header", require_arg, 'H' },  // Print shorter version without the values
    { NULL, 0, '\0' }
};

/* a structure for handling the order command-line parameters come in */
struct arg_handler {
    void (*func)(h5_file_t *, char *);
    char *obj;
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
static struct arg_handler* function_assign(int argc, const char *argv[])
{
    struct arg_handler   *hand = NULL;
 
    int                  i, option;

    /* this will be plenty big enough to hold the info */
    hand = (arg_handler*)calloc((size_t)argc, sizeof(struct arg_handler));

    /* set options according to the command line */
    while ((option = get_option(argc, argv, s_opts, l_opts)) != EOF)
    {
       switch ((char)option)
       {
          case 'h': // Print help page
            print_help();
            exit(1);
          case 'A': // Print file attributes
            display_all = 0;

            for (i = 0; i < argc; i++)
            {
               if (!hand[i].func)
               {
                  hand[i].func = print_file_attributes;
                  hand[i].obj = NULL; // inserting garabage value that we won't use. (For function interface compatibility)
                  break;
               }
            }
            break;
          case 'a': // Print step attributes & values for time step n
            display_all = 0;

            for (i = 0; i < argc; i++)
            {
               if (!hand[i].func)
               {
                  hand[i].func = print_step_attributes;
                  hand[i].obj = strdup(opt_arg);
                  break;
               }
            }
            break;
          case 'd': // Print data sets names & values for time step n
            display_all = 0;

            for (i = 0; i < argc; i++)
            {
               if (!hand[i].func)
               {
                  hand[i].func = print_dataset;
                  hand[i].obj = strdup(opt_arg);
                  break;
               }
            }
            break;
          case 'n': // Print number of steps
            display_all = 0;

            for (i = 0; i < argc; i++)
            {
               if (!hand[i].func)
               {
                  hand[i].func = print_nstep;
                  hand[i].obj = NULL; // inserting garabage value that we won't use. (For function interface compatibility)
                  break;
               }
            }
          break;
          case 'H': // Print shorter version without the values
            print_header = true;
            break;
          default:
            print_help();
            exit(1);
       }
    }
    return hand;
}

/* For printing help page */
static void print_help()
{
   fflush(stdout);
   fprintf(stdout, "\nusage: h5pAttrib [OPTIONS] file\n");
   fprintf(stdout, "\n");
   fprintf(stdout, "  OPTIONS\n");
   fprintf(stdout, "   -h, --help           Print help page\n");
   fprintf(stdout, "   -n, --nstep          Print number of steps\n");
   fprintf(stdout, "   -A, --fileA          Print file attributes\n");
   fprintf(stdout, "   -a n, --stepA n      Print step attributes & values for time step n\n");
   fprintf(stdout, "   -d n, --dataset n    Print data sets names & values for time step n\n");
   fprintf(stdout, "   -H, --header         Print shorter version without the values\n");
   fprintf(stdout, "\n");
   fprintf(stdout, "  Examples:\n");
   fprintf(stdout, "\n");
   fprintf(stdout, "  1) Show file attribute names & values of sample.h5part\n");
   fprintf(stdout, "\n");
   fprintf(stdout, "        h5pAttrib -A sample.h5part\n");
   fprintf(stdout, "\t\t\tOR\n");
   fprintf(stdout, "        h5pAttrib --fileA sample.h5part\n");
   fprintf(stdout, "\n");
   fprintf(stdout, "  2) Show step attribute names for time step 5 of sample.h5part\n");
   fprintf(stdout, "\n");
   fprintf(stdout, "        h5pAttrib -a 5 -H sample.h5part\n");
   fprintf(stdout, "\t\t\tOR\n");
   fprintf(stdout, "        h5pAttrib --stepA 5 -H sample.h5part\n");
   fprintf(stdout, "\n");
}

/* For priting everything (default option when no flags provided.) */
static void print_all(h5_file_t* file)
{
   int nt;

   char file_attrib[MAX_LEN];
   h5_int64_t type;
   h5_int64_t num_elem;
   h5_int64_t num_attrib;
   char step_attrib[MAX_LEN];
   h5_int64_t i, j, k;
   h5_int64_t count;
   void* value = NULL;
   h5_int64_t timestep_ctr;

   char data_name[MAX_LEN];
   h5_int64_t num_dataset;
   h5_int64_t nparticles;

   fprintf(stdout, "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
   fprintf(stdout, "\n");
   nt=H5GetNumSteps(file);
   fprintf(stdout, "There are total %d number of timesteps.\n",nt);

   fprintf(stdout, "\nDump result for \"%s\"...\n\n", global_fname);
   fprintf(stdout, "FILE_ATTRIBUTES:\n");


   if(print_header)
   {
      num_attrib = H5GetNumFileAttribs(file);

      for(i=0; i<num_attrib; i++)
      {
         H5GetFileAttribInfo(file, i, file_attrib, MAX_LEN, &type, &num_elem);
         fprintf(stdout, "\tAttribute #%lld = %s\n",(long long)i,file_attrib);
         fprintf(stdout, "\tThere are %lld elements in the attribute\n",(long long)num_elem);

         if (type == H5_INT64_T)
         {
            fprintf(stdout, "\tAttribute Type: H5_INT64_T\n");
         }
         else if (type == H5_INT32_T)
         {
            fprintf(stdout, "\tAttribute Type: H5_INT32_T\n");
         }
         else if (type == H5_FLOAT64_T)
         {
            fprintf(stdout, "\tAttribute Type: H5_FLOAT64_T\n");
         }
         else if (type == H5_FLOAT32_T)
         {
            fprintf(stdout, "\tAttribute Type: H5_FLOAT32_T\n");
         }
         else if (type == H5_STRING_T)
         {
            fprintf(stdout, "\tAttribute Type: H5_STRING_T\n");
         }
         else
         {
             fprintf(stdout, "\tAttribute Type is UNKNOWN.\n");
         }
      }
      if(i==0)
      {
         fprintf(stdout, "\tNO FILE ATTRIBUTES...\n");
      }
   }
   else
   {
      num_attrib = H5GetNumFileAttribs(file);

      for(i=0; i<num_attrib; i++)
      {
         H5GetFileAttribInfo(file, i, file_attrib, MAX_LEN, &type, &num_elem);
         fprintf(stdout, "\tAttribute #%lld = %s\n",(long long)i,file_attrib);
         fprintf(stdout, "\tThere are %lld elements in the attribute\n",(long long)num_elem);

         if (type == H5T_NATIVE_INT64)
         {
            value = (h5_int64_t*)malloc(sizeof(h5_int64_t)*num_elem);
            H5PartReadFileAttrib(file, file_attrib, (h5_int64_t*)value);

            fprintf(stdout, "\tAttribute Type: H5T_NATIVE_INT64\n");

            fprintf(stdout, "\tPrinting %lld element value(s):\n", (long long)num_elem);

            k=0;
            count=1;
            fprintf(stdout, "\tCOUNT[0]: ");
            for(j=0; j<num_elem; j++)
            {
               fprintf(stdout, "%lld  ",((long long*)value)[j]);
               k++;
               if(k==5)
               {
                  fprintf(stdout, "\n");
                  fprintf(stdout, "\tCOUNT[%lld]: ", (long long)count);
                  k=0;
               }
               count++;
            }
            fprintf(stdout, "\n");

            free(value);
         }
         else if (type == H5T_NATIVE_CHAR)
         {
            value = (char*)malloc(sizeof(char)*num_elem);
            H5PartReadFileAttrib(file, file_attrib, (char*)value);

            fprintf(stdout, "\tAttribute Type: H5T_NATIVE_CHAR\n");

            fprintf(stdout, "\tPrinting string of length %lld:\n", (long long)num_elem);

            k=0;
            for(j=0; j<num_elem; j++)
            {
               fprintf(stdout, "%c",((char*)value)[j]);
               k++;
               if(k==100)
               {
                  fprintf(stdout, "\n");
                  k=0;
               }
            }
            fprintf(stdout, "\n");

            free(value);
         }
         else if (type == H5T_NATIVE_DOUBLE)
         {
            value = (double*)malloc(sizeof(double)*num_elem);
            H5PartReadFileAttrib(file, file_attrib, (double*)value);

            fprintf(stdout, "\tAttribute Type: H5T_NATIVE_DOUBLE\n");

            fprintf(stdout, "\tPrinting %lld element value(s):\n", (long long)num_elem);

            k=0;
            count=1;
            fprintf(stdout, "\tCOUNT[0]: ");
            for(j=0; j<num_elem; j++)
            {
               fprintf(stdout, "%lf  ", ((double*)value)[j]);
               k++;
               if(k==5)
               {
                  fprintf(stdout, "\n");
                  fprintf(stdout, "\tCOUNT[%lld]: ", (long long)count);
                  k=0;
               }
               count++;
            }

            fprintf(stdout, "\n");

            free(value);
         }
         else
         {
             fprintf(stdout, "\tAttribute Type is UNKNOWN.\n");
         }
         fprintf(stdout, "\n");
      }
      if(i==0)
      {
         fprintf(stdout, "\tNO FILE ATTRIBUTES...\n");
      }
   }
   fprintf(stdout, "\n");

   for(timestep_ctr=0; timestep_ctr<nt; timestep_ctr++)
   {
       H5PartSetStep(file,timestep_ctr);

       fprintf(stdout, "Timestep #%lld:\n", (long long)timestep_ctr);
       num_attrib = H5PartGetNumStepAttribs(file);

       if(print_header)
       {
          for(i=0; i<num_attrib; i++)
          {
             H5PartGetStepAttribInfo(file, i, step_attrib, MAX_LEN, &type, &num_elem);
             fprintf(stdout, "\tAttribute #%lld = %s\n",(long long)i,step_attrib);
             fprintf(stdout, "\tNumber of elements in the attribute: %lld\n",(long long)num_elem);

             if (type == H5T_NATIVE_INT64)
             {
                fprintf(stdout, "\tAttribute Type: H5T_NATIVE_INT64\n");
             }
             else if (type == H5T_NATIVE_CHAR)
             {
                fprintf(stdout, "\tAttribute Type: H5T_NATIVE_CHAR\n");
             }
             else if (type == H5T_NATIVE_DOUBLE)
             {
                fprintf(stdout, "\tAttribute Type: H5T_NATIVE_DOUBLE\n");
             }
             else
             {
                fprintf(stdout, "\tAttribute Type is UNKNOWN.\n");
             }
             fprintf(stdout, "\n");
          }
          if(i==0)
          {
             fprintf(stdout, "\tNO STEP ATTRIBUTES...\n");
          }
       }
       else
       {
          for(i=0; i<num_attrib; i++)
          {
             H5PartGetStepAttribInfo(file, i, step_attrib, MAX_LEN, &type, &num_elem);
             fprintf(stdout, "\tAttribute #%lld = %s\n",(long long)i,step_attrib);
             fprintf(stdout, "\tNumber of elements in the attribute: %lld\n",(long long)num_elem);

             if (type == H5T_NATIVE_INT64)
             {
                value = (long long int*)malloc(sizeof(h5_int64_t)*num_elem);
                H5PartReadStepAttrib(file, step_attrib, (h5_int64_t*)value);

                fprintf(stdout, "\tAttribute Type: H5T_NATIVE_INT64\n");

                fprintf(stdout, "\tPrinting %lld element value(s):\n", (long long)num_elem);

                k=0;
                count=1;
                fprintf(stdout, "\tCOUNT[0]: ");
                for(j=0; j<num_elem; j++)
                {
                   fprintf(stdout, "%lld  ",((long long*)value)[j]);
                   k++;
                   if(k==5)
                   {
                      fprintf(stdout, "\n");
                      fprintf(stdout, "\tCOUNT[%lld]: ", (long long)count);
                      k=0;
                   }
                   count++;
                }
                fprintf(stdout, "\n");

                free(value);
             }
             else if (type == H5T_NATIVE_CHAR)
             {
                value = (char*)malloc(sizeof(char)*num_elem);
                H5PartReadStepAttrib(file, step_attrib, (char*)value);

                fprintf(stdout, "\tAttribute Type: H5T_NATIVE_CHAR\n");

                fprintf(stdout, "\tPrinting string of length %lld:\n", (long long)num_elem);

                k=0;
                for(j=0; j<num_elem; j++)
                {
                   fprintf(stdout, "%c",((char*)value)[j]);
                   k++;
                   if(k==100)
                   {
                      fprintf(stdout, "\n");
                      k=0;
                   }
                }
                fprintf(stdout, "\n");

                free(value);
             }
             else if (type == H5T_NATIVE_DOUBLE)
             {
                value = (double*)malloc(sizeof(double)*num_elem);
                H5PartReadStepAttrib(file, step_attrib, (double*)value);

                fprintf(stdout, "\tAttribute Type: H5T_NATIVE_DOUBLE\n");

                fprintf(stdout, "\tPrinting %lld element value(s):\n", (long long)num_elem);

                k=0;
                count=1;
                fprintf(stdout, "\tCOUNT[0]: ");
                for(j=0; j<num_elem; j++)
                {
                   fprintf(stdout, "%lf  ", ((double*)value)[j]);
                   k++;
                   if(k==5)
                   {
                      fprintf(stdout, "\n");
                      fprintf(stdout, "\tCOUNT[%lld]: ", (long long)count);
                      k=0;
                   }
                   count++;
                }

                fprintf(stdout, "\n");

                free(value);
             }
             else
             {
                fprintf(stdout, "\tAttribute Type is UNKNOWN.\n");
             }
             fprintf(stdout, "\n");
          }
          if(i==0)
          {
             fprintf(stdout, "\tNO STEP ATTRIBUTES...\n\n");
          }
       }

      num_dataset = H5PartGetNumDatasets(file);
      if(print_header)
      {
         for(i=0; i<num_dataset; i++)
         {
            H5PartGetDatasetInfo(file, i, data_name, MAX_LEN, &type, &nparticles);
            fprintf(stdout, "\tDataset Name #%lld = %s\n",(long long)i,data_name);
            fprintf(stdout, "\tNumber of elements in the dataset: %lld\n",(long long)nparticles);

            if (type == H5T_NATIVE_INT64)
            {
               fprintf(stdout, "\tDataset Type: H5T_NATIVE_INT64\n");
            }
            else if (type == H5T_NATIVE_DOUBLE)
            {
               fprintf(stdout, "\tDataset Type: H5T_NATIVE_DOUBLE\n");
            }
            else
            {
               fprintf(stdout, "\tDataset Type: UNKNOWN.\n");
            }
            fprintf(stdout, "\n");
         }

         if(i==0)
         {
             fprintf(stdout, "\tNO DATASETS...\n\n");
         }
      }
      else
      {
         for(i=0; i<num_dataset; i++)
         {
            H5PartGetDatasetInfo(file, i, data_name, MAX_LEN, &type, &nparticles);
            fprintf(stdout, "\tDataset Name #%lld = %s\n",(long long)i,data_name);

            if (type == H5T_NATIVE_INT64)
            {
               value = (h5_int64_t*)malloc(sizeof(h5_int64_t)*nparticles);
               H5PartReadDataInt64(file, data_name, (h5_int64_t*)value);
               fprintf(stdout, "\tDataset Type: H5T_NATIVE_INT64\n");
               fprintf(stdout, "\tPrinting %lld element value(s):\n", (long long)nparticles);

               count=1;
               k=0;
               fprintf(stdout, "\tCOUNT[0]: ");
               for(j=0; j<nparticles; j++)
               {
                  fprintf(stdout, "%lld  ", ((long long*)value)[j]);
                  k++;
                  if(k==5)
                  {
                     fprintf(stdout, "\n");
                     fprintf(stdout, "\tCOUNT[%lld]: ", (long long)count);
                     k=0;
                  }
                  count++;
               }

               fprintf(stdout, "\n");

               free(value);
            }
            else if (type == H5T_NATIVE_DOUBLE)
            {
               value = (double*)malloc(sizeof(double)*nparticles);
               H5PartReadDataFloat64(file, data_name, (double*)value);
               fprintf(stdout, "\tDataset Type: H5T_NATIVE_DOUBLE\n");
               fprintf(stdout, "\tPrinting %lld element value(s):\n", (long long)nparticles);

               k=0;
               count=1;
               fprintf(stdout, "\tCOUNT[0]: ");
               for(j=0; j<nparticles; j++)
               {
                  fprintf(stdout, "%lf  ", ((double*)value)[j]);
                  k++;
                  if(k==5)
                  {
                     fprintf(stdout, "\n");
                     fprintf(stdout, "\tCOUNT[%lld]: ", (long long)count);
                     k=0;
                  }
                  count++;
               }

               fprintf(stdout, "\n");

               free(value);
            }
            else
            {
               fprintf(stdout, "\tDataset Type: UNKNOWN.\n");
               fprintf(stdout, "\n");
            }
            fprintf(stdout, "\n");
         }

         if(i==0)
         {
             fprintf(stdout, "\tNO DATASETS...\n\n");
         }
      }
   }
    fprintf(stdout, "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n");
}

/* For printing number of timesteps in the file */
static void print_nstep(h5_file_t* file, char * garbage)
{
    int nt;

    fprintf(stdout, "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    fprintf(stdout, "\nPrinting number of timesteps for: %s ...\n", global_fname);
    fprintf(stdout, "\n");
    nt=H5PartGetNumSteps(file);
    fprintf(stdout, "There are total %d number of timesteps.\n",nt);
    fprintf(stdout, "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n");
}

/* For printing file attributes */
static void print_file_attributes(h5_file_t* file, char * garbage)
{
    char file_attrib[MAX_LEN];
    h5_int64_t type;
    h5_int64_t num_elem;
    h5_int64_t num_attrib;
    h5_int64_t i, j, k;
    h5_int64_t count;
    void* value = NULL;

    if(print_header)
    {
       fprintf(stdout, "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
       fprintf(stdout, "\nPrinting file attributes for: %s ...\n", global_fname);
       fprintf(stdout, "\n");
       num_attrib = H5PartGetNumFileAttribs(file);
       fprintf(stdout, "The number of file attributes for file %s is %lld ...\n", global_fname, (long long)num_attrib);

       for(i=0; i<num_attrib; i++)
       {
          H5PartGetFileAttribInfo(file, i, file_attrib, MAX_LEN, &type, &num_elem);
          fprintf(stdout, "Attribute #%lld = %s\n",(long long)i,file_attrib);

          if (type == H5T_NATIVE_INT64)
          {
             fprintf(stdout, "Attribute Type: H5T_NATIVE_INT64\n");
          }
          else if (type == H5T_NATIVE_CHAR)
          {
             fprintf(stdout, "Attribute Type: H5T_NATIVE_CHAR\n");
          }
          else if (type == H5T_NATIVE_DOUBLE)
          {
             fprintf(stdout, "Attribute Type: H5T_NATIVE_DOUBLE\n");
          }
          else
          {
              fprintf(stdout, "Attribute Type is UNKNOWN.\n");
          }
       }
       if(i==0)
       {
           fprintf(stdout, "There are no file attributes.\n");
       }
       fprintf(stdout, "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n");
    }
    else
    {
       fprintf(stdout, "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
       fprintf(stdout, "\nPrinting file attributes for: %s ...\n", global_fname);
       fprintf(stdout, "\n");
       num_attrib = H5PartGetNumFileAttribs(file);
       fprintf(stdout, "The number of file attributes for file %s is %lld ...\n", global_fname, (long long)num_attrib);

       for(i=0; i<num_attrib; i++)
       {
          H5PartGetFileAttribInfo(file, i, file_attrib, MAX_LEN, &type, &num_elem);
          fprintf(stdout, "Attribute #%lld = %s\n",(long long)i,file_attrib);
          fprintf(stdout, "There are %lld elements in the attribute\n",(long long)num_elem);

          if (type == H5T_NATIVE_INT64)
          {
             value = (h5_int64_t*)malloc(sizeof(h5_int64_t)*num_elem);
             H5PartReadFileAttrib(file, file_attrib, (h5_int64_t*)value);

             fprintf(stdout, "Attribute Type: H5T_NATIVE_INT64\n");

             fprintf(stdout, "Printing %lld element value(s):\n", (long long)num_elem);

             k=0;
             count=1;
             fprintf(stdout, "COUNT[0]: ");
             for(j=0; j<num_elem; j++)
             {
                fprintf(stdout, "%lld  ",((long long*)value)[j]);
                k++;
                if(k==5)
                {
                   fprintf(stdout, "\n");
                   fprintf(stdout, "COUNT[%lld]: ", (long long)count);
                   k=0;
                }
                count++;
             }
             fprintf(stdout, "\n");

             free(value);
          }
          else if (type == H5T_NATIVE_CHAR)
          {
             value = (char*)malloc(sizeof(char)*num_elem);
             H5PartReadFileAttrib(file, file_attrib, (char*)value);

             fprintf(stdout, "Attribute Type: H5T_NATIVE_CHAR\n");

             fprintf(stdout, "Printing string of length %lld:\n", (long long)num_elem);

             k=0;
             for(j=0; j<num_elem; j++)
             {
                fprintf(stdout, "%c",((char*)value)[j]);
                k++;
                if(k==100)
                {
                   fprintf(stdout, "\n");
                   k=0;
                }
             }
             fprintf(stdout, "\n");

             free(value);
          }
          else if (type == H5T_NATIVE_DOUBLE)
          {
             value = (double*)malloc(sizeof(double)*num_elem);
             H5PartReadFileAttrib(file, file_attrib, (double*)value);

             fprintf(stdout, "Attribute Type: H5T_NATIVE_DOUBLE\n");

             fprintf(stdout, "Printing %lld element value(s):\n", (long long)num_elem);

             k=0;
             count=1;
             fprintf(stdout, "COUNT[0]: ");
             for(j=0; j<num_elem; j++)
             {
                fprintf(stdout, "%lf  ", ((double*)value)[j]);
                k++;
                if(k==5)
                {
                   fprintf(stdout, "\n");
                   fprintf(stdout, "COUNT[%lld]: ", (long long)count);
                   k=0;
                }
                count++;
             }

             fprintf(stdout, "\n");

             free(value);
          }
          else
          {
              fprintf(stdout, "Attribute Type is UNKNOWN.\n");
          }
          fprintf(stdout, "\n");
       }
       if(i==0)
       {
           fprintf(stdout, "There are no file attributes.\n");
       }
       fprintf(stdout, "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n");
    }
}

/* For printing step attribute names [& values] */
static void print_step_attributes(h5_file_t* file, char *attr)
{
    char step_attrib[MAX_LEN];
    h5_int64_t type;
    h5_int64_t num_elem;
    h5_int64_t num_attrib;
    h5_int64_t i, j, k;
    h5_int64_t count;
    void* value = NULL;

    if(print_header)
    {
       fprintf(stdout, "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
       fprintf(stdout, "\nPrinting step attributes for: %s ...\n", global_fname);
       fprintf(stdout, "\n");
       H5PartSetStep(file,atoi(attr));
       num_attrib = H5PartGetNumStepAttribs(file);
       fprintf(stdout, "The number of step attributes for timestep #%d is %lld ...\n\n", atoi(attr), (long long)num_attrib);
       for(i=0; i<num_attrib; i++)
       {
          H5PartGetStepAttribInfo(file, i, step_attrib, MAX_LEN, &type, &num_elem);
          fprintf(stdout, "Attribute #%lld = %s\n",(long long)i,step_attrib);
       }
       if(i==0)
       {
           fprintf(stdout, "There are no step attributes for timestep #%d.\n", atoi(attr));
       }
       fprintf(stdout, "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n");

    }
    else
    {
       fprintf(stdout, "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
       fprintf(stdout, "\nPrinting step attributes for: %s ...\n", global_fname);
       fprintf(stdout, "\n");
       H5PartSetStep(file,atoi(attr));
       num_attrib = H5PartGetNumStepAttribs(file);
       fprintf(stdout, "The number of step attributes for timestep #%d is %lld ...\n\n", atoi(attr), (long long)num_attrib);
       for(i=0; i<num_attrib; i++)
       {
          H5PartGetStepAttribInfo(file, i, step_attrib, MAX_LEN, &type, &num_elem);
          fprintf(stdout, "Attribute #%lld = %s\n",(long long)i,step_attrib);
          fprintf(stdout, "There are %lld elements in the attribute\n",(long long)num_elem);

          if (type == H5T_NATIVE_INT64)
          {
             value = (h5_int64_t*)malloc(sizeof(h5_int64_t)*num_elem);
             H5PartReadStepAttrib(file, step_attrib, (h5_int64_t*)value);

             fprintf(stdout, "Attribute Type is H5T_NATIVE_INT64\n");

             fprintf(stdout, "Printing %lld element value(s):\n", (long long)num_elem);

             k=0;
             count=1;
             fprintf(stdout, "COUNT[0]: ");
             for(j=0; j<num_elem; j++)
             {
                fprintf(stdout, "%lld  ",((long long*)value)[j]);
                k++;
                if(k==5)
                {
                   fprintf(stdout, "\n");
                   fprintf(stdout, "COUNT[%lld]: ", (long long)count);
                   k=0;
                }
                count++;
             }
             fprintf(stdout, "\n");

             free(value);
          }
          else if (type == H5T_NATIVE_CHAR)
          {
             value = (char*)malloc(sizeof(char)*num_elem);
             H5PartReadStepAttrib(file, step_attrib, (char*)value);

             fprintf(stdout, "Attribute Type is H5T_NATIVE_CHAR\n");

             fprintf(stdout, "Printing string of length %lld:\n", (long long)num_elem);

             k=0;
             for(j=0; j<num_elem; j++)
             {
                fprintf(stdout, "%c",((char*)value)[j]);
                k++;
                if(k==100)
                {
                   fprintf(stdout, "\n");
                   k=0;
                }
             }
             fprintf(stdout, "\n");

             free(value);
          }
          else if (type == H5T_NATIVE_DOUBLE)
          {
             value = (double*)malloc(sizeof(double)*num_elem);
             H5PartReadStepAttrib(file, step_attrib, (double*)value);

             fprintf(stdout, "Attribute Type is H5T_NATIVE_DOUBLE\n");

             fprintf(stdout, "Printing %lld element value(s):\n", (long long)num_elem);

             k=0;
             count=1;
             fprintf(stdout, "COUNT[0]: ");
             for(j=0; j<num_elem; j++)
             {
                fprintf(stdout, "%lf  ", ((double*)value)[j]);
                k++;
                if(k==5)
                {
                   fprintf(stdout, "\n");
                   fprintf(stdout, "COUNT[%lld]: ", (long long)count);
                   k=0;
                }
                count++;
             }

             fprintf(stdout, "\n");

             free(value);
          }
          else
          {
              fprintf(stdout, "Attribute Type is UNKNOWN.\n");
          }
          fprintf(stdout, "\n");
       }
       if(i==0)
       {
           fprintf(stdout, "There are no step attributes for timestep #%d.\n", atoi(attr));
       }
       fprintf(stdout, "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n");
    }
}

/* For printing dataset names [& values] */
static void print_dataset(h5_file_t* file, char *attr)
{
    char data_name[MAX_LEN];
    int i, j, k;
    long count;
    h5_int64_t type;
    int num_dataset;
    h5_int64_t nparticles;
    void* value = NULL;


    if(print_header)
    {
       fprintf(stdout, "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
       fprintf(stdout, "\nPrinting names of datasets for: %s ...\n", global_fname);
       fprintf(stdout, "\n");
       H5PartSetStep(file,atoi(attr));
       num_dataset = H5PartGetNumDatasets(file);
       fprintf(stdout, "The number of datasets for timestep #%d is %d ...\n\n", atoi(attr), num_dataset);

       for(i=0; i<num_dataset; i++)
       {
          H5PartGetDatasetInfo(file, i, data_name, MAX_LEN, &type, &nparticles);
          fprintf(stdout, "Dataset Name #%d = %s\n",i,data_name);

          if (type == H5T_NATIVE_INT64)
          {
             fprintf(stdout, "Dataset Type is H5T_NATIVE_INT64\n");
             fprintf(stdout, "Number of elements: %lld\n", (long long)nparticles);
          }
          else if (type == H5T_NATIVE_DOUBLE)
          {
             fprintf(stdout, "Dataset Type is H5T_NATIVE_DOUBLE\n");
             fprintf(stdout, "Number of elements: %lld\n", (long long)nparticles);
          }
          else
          {
             fprintf(stdout, "Dataset Type is UNKNOWN.\n");
          }
          fprintf(stdout, "\n");
       }

       if(i==0)
       {
           fprintf(stdout, "There are no datasets for timestep #%d.\n", atoi(attr));
       }
       fprintf(stdout, "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n");

    }
    else
    {
       fprintf(stdout, "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
       fprintf(stdout, "\nPrinting names of datasets for: %s ...\n", global_fname);
       fprintf(stdout, "\n");
       H5SetStep(file,atoi(attr));
       num_dataset = H5PartGetNumDatasets(file);
       fprintf(stdout, "The number of datasets for timestep #%d is %d ...\n\n", atoi(attr), num_dataset);

       for(i=0; i<num_dataset; i++)
       {
          H5PartGetDatasetInfo(file, i, data_name, MAX_LEN, &type, &nparticles);
          fprintf(stdout, "Dataset Name #%d = %s\n",i,data_name);

          if (type == H5T_NATIVE_INT64)
          {
             value = (h5_int64_t*)malloc(sizeof(h5_int64_t)*nparticles);
             H5PartReadDataInt64(file, data_name, (h5_int64_t*)value);
             fprintf(stdout, "Dataset Type is H5T_NATIVE_INT64\n");
             fprintf(stdout, "Printing %lld element value(s):\n", (long long)nparticles);

             count=1;
             k=0;
             fprintf(stdout, "COUNT[0]: ");
             for(j=0; j<nparticles; j++)
             {
                fprintf(stdout, "%lld  ", ((long long*)value)[j]);
                k++;
                if(k==5)
                {
                   fprintf(stdout, "\n");
                   fprintf(stdout, "COUNT[%ld]: ", count);
                   k=0;
                }
                count++;
             }

             fprintf(stdout, "\n");

             free(value);
          }
          else if (type == H5T_NATIVE_DOUBLE)
          {
             value = (double*)malloc(sizeof(double)*nparticles);
             H5PartReadDataFloat64(file, data_name, (double*)value);
             fprintf(stdout, "Dataset Type is H5T_NATIVE_DOUBLE\n");
             fprintf(stdout, "Printing %lld element value(s):\n", (long long)nparticles);

             k=0;
             count=1;
             fprintf(stdout, "COUNT[0]: ");
             for(j=0; j<nparticles; j++)
             {
                fprintf(stdout, "%lf  ", ((double*)value)[j]);
                k++;
                if(k==5)
                {
                   fprintf(stdout, "\n");
                   fprintf(stdout, "COUNT[%ld]: ", count);
                   k=0;
                }
                count++;
             }

             fprintf(stdout, "\n");

             free(value);
          }
          else
          {
             fprintf(stdout, "Dataset Type is UNKNOWN.\n");
             fprintf(stdout, "\n");
          }
          fprintf(stdout, "\n");
       }

       if(i==0)
       {
           fprintf(stdout, "There are no datasets for timestep #%d.\n", atoi(attr));
       }

       fprintf(stdout, "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n");
    }
}

/* Frees argument handlers */
static void free_handler(struct arg_handler *hand, int len)
{
    int i;

    for (i = 0; i < len; i++) 
    {
        free(hand[i].obj);
    }

    free(hand);
}

int main(int argc, const char *argv[])
{
   /* Numerous variables */
   struct arg_handler   *hand = NULL;
   int                   i;
   h5_file_t           *h5file = NULL;
   const char         *fname = NULL;

   //h5pAttrib_function_table = &major_function_table;
   /* Take care of the command line options */
   hand = function_assign(argc, argv);

    if (argc <= opt_ind) 
    {
        fprintf(stdout, "missing file name\n");
        print_help();
        exit(1);
    }


   /* Check for conflicting options */
   /* There are none. If -H is appended to non-compatible flags, just ignore. */


   /* Process accordingly */
   fname = argv[opt_ind];
   global_fname = fname; // To use in funtions
   h5file = H5OpenFile(fname,H5_O_RDONLY,0);

   if ( h5file == NULL )
   {
        fprintf(stdout, "unable to open file %s\n", fname);
        print_help();
        exit(1);
   }

   if (display_all)
   {
      print_all(h5file);
   }
   else 
   {
      for (i = 0; i < argc; i++)
      {
         if (hand[i].func)
         {
            hand[i].func(h5file, hand[i].obj);
         }
      }
   }

   free_handler(hand, argc);
   H5CloseFile(h5file);
   fprintf(stdout,"done\n");
   return 0;
}
