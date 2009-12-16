/* Command line interface for the fldiagnostics module
 *
 * James Maddison
 *
 * Command line argument handling based on ADmain.cpp
 */

#include <assert.h>
#include <getopt.h>
#include <iostream>
#include <map>
#include <queue>
#include <stdio.h>
#include <string>

#include "confdefs.h"
#include "spud.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

/*
#define fldiag_diff_fc F77_FUNC(fldiag_diff, FLDIAG_DIFF)
extern "C"{
  void fldiag_diff_fc(const char* input_name, const int* input_name_len,
                      const char* output_name, const int* output_name_len,
                      const char* infield_name, const int* infield_name_len,
                      const char* outfield_name, const int* outfield_name_len,
                      const char* diff_type, const int* diff_type_len,
                      const char* coords_name, const int* coords_name_len);
}
*/

#define fldiag_add_diag_fc F77_FUNC(fldiag_add_diag, FLDIAG_ADD_DIAG)
extern "C"{
  void fldiag_add_diag_fc(const char* input_name, const int* input_name_len,
                          const char* output_name, const int* output_name_len,
                          const char* outfield_name,
                            const int* outfield_name_len,
                          const char* meshfield_name,
                            const int* meshfield_name_len,
                          const char* state_name,
                            const int* state_name_len,
                          const int* outfield_rank);
}

/*
#define fldiag_add_presc_fc F77_FUNC(fldiag_add_presc, FLDIAG_ADD_PRESC)
extern "C"{
  void fldiag_add_presc_fc(const char* input_name, const int* input_name_len,
                           const char* output_name, const int* output_name_len,
                           const char* outfield_name,
                             const int* outfield_name_len,
                           const char* meshfield_name,
                             const int* meshfield_name_len,
                           const char* state_name,
                             const int* state_name_len);
}
*/

void Help();
bool ParseArgs(int argc, char** argv);
bool HaveOption(const std::string& opt);
const char* GetOption(const std::string& opt);
const char* GetOption(const std::string& opt, int& stat);
void FatalError(const char* message);
int main(int argc, char** argv);

void Help(){
  std::cout << "Usage: fldiagnostics ACTION [OPTIONS] INPUT OUTPUT [FIRST] [LAST]\n"
            << "\n"
            << "Offline diagnostics tools.\n"
            << "\n"
            << "If FIRST is supplied, treats INPUT and OUTPUT as project names, and processes\n"
            << "the specified range of project files.\n"
            << "\n"
            << "Actions:\n"
            << "\n"
            << "add[-diag]  Add a diagnostics field. Options:\n"
            << "              -m NAME   Field from which to extract a mesh to use with the\n"
            << "                        diagnostic field (default \"Velocity\")\n"
            << "              -o NAME   Diagnostic field name\n"
            << "              -r RANK   Specify the rank of the diagnostic field (will try all\n"
            << "                        ranks if not supplied, but will also suppress useful\n"
            << "                        error messages)\n"
            << "              -s STATE  Name of the state from which to read options in the \n"
            << "                        options file (defaults to the first state)\n"
            << "              -x FLML   Model options file (not always required)"
/*          << "\n"
            << "add-presc   Add a prescribed field, as specified from the supplied options file.\n"
            << "            Options:\n"
            << "              -m NAME   Field from which to extract a mesh to use with the\n"
            << "                        prescribed field (default \"Velocity\")\n"
            << "              -o NAME   Prescribed field name\n"
            << "              -s STATE  Name of the state from which to read options in the \n"
            << "                        options file (defaults to the first state)\n"
            << "              -x FLML   Model options file"
            << "\n"
            << "diff  Take the derivative of the supplied field. Options:\n"
            << "        -d TYPE  Derivative type. Supported types:\n"
            << "                   grad  Gradient operator\n"
            << "        -i NAME  Input field name\n"
            << "        -o NAME  Output field name\n"
            << "        -p NAME  Set the name of the coordinates field
  (default \"Coordinate\")\n"
*/
            << std::endl;

  return;
}

std::map< std::string, std::string > fldiagnostics_opts;

bool ParseArgs(int argc, char** argv){
  if(argc < 2){
    return false;
  }
  fldiagnostics_opts["action"] = argv[1];

  opterr = 0;
  optind = 2;
  if((strlen(argv[1]) == 8 and strncmp(argv[1], "add-diag", 8) == 0) or (strlen(argv[1]) == 3 and strncmp(argv[1], "add", 3) == 0)){
    fldiagnostics_opts["meshfield_name"] = "Velocity";

    char opt;
    while(true){
      opt = getopt(argc, argv, "m:o:r:s:x:");
      if(opt == '?'){
        return false;
      }
      else if(opt == -1){
        break;
      }
      switch(opt){
        case 'm':
          if(optarg == NULL){
            return false;
          }
          fldiagnostics_opts["meshfield_name"] = optarg;
          break;

        case 'o':
          if(optarg == NULL){
            return false;
          }
          fldiagnostics_opts["outfield_name"] = optarg;
          break;

        case 'r':
          int test;
          try{
            test = atoi(optarg);
          }
          catch(...){
            return false;
          }
          if(test < 0 or test > 2){
            return false;
          }
          fldiagnostics_opts["outfield_rank"] = optarg;
          break;

        case 's':
          if(optarg == NULL){
            return false;
          }
          fldiagnostics_opts["state_name"] = optarg;
          break;

        case 'x':
          if(optarg == NULL){
            return false;
          }
          fldiagnostics_opts["options_file"] = optarg;
          break;

        default:
          break;
      }
    }
    if(!HaveOption("outfield_name")){
      return false;
    }
  }
/*
  else if(strlen(argv[1]) == 9 and strncmp(argv[1], "add-presc", 9) == 0){
    fldiagnostics_opts["meshfield_name"] = "Velocity";

    char opt;
    while(true){
      opt = getopt(argc, argv, "m:o:r:s:x:");
      if(opt == '?'){
        return false;
      }
      else if(opt == -1){
        break;
      }
      switch(opt){
        case 'm':
          if(optarg == NULL){
            return false;
          }
          fldiagnostics_opts["meshfield_name"] = optarg;
          break;

        case 'o':
          if(optarg == NULL){
            return false;
          }
          fldiagnostics_opts["outfield_name"] = optarg;
          break;

        case 's':
          if(optarg == NULL){
            return false;
          }
          fldiagnostics_opts["state_name"] = optarg;
          break;

        case 'x':
          if(optarg == NULL){
            return false;
          }
          fldiagnostics_opts["options_file"] = optarg;
          break;

        default:
          break;
      }
    }
    if(!HaveOption("outfield_name") or !HaveOption("options_file")){
      return false;
    }
  }
*/
  else{
    return false;
  }

  if(optind + 1 >= argc or optind + 5 <= argc){
    return false;
  }

  fldiagnostics_opts["input_name"] = argv[optind];
  fldiagnostics_opts["output_name"] = argv[optind + 1];
  if(optind + 2 < argc){
    try{
      (unsigned int)atoi(argv[optind  + 2]);
    }
    catch(...){
      return false;
    }
    fldiagnostics_opts["first_id"] = argv[optind  + 2];
    if(optind + 3 < argc){
      try{
        (unsigned int)atoi(argv[optind + 3]);
      }
      catch(...){
        return false;
      }
      fldiagnostics_opts["last_id"] = argv[optind + 3];
    }
    else{
      fldiagnostics_opts["last_id"] = argv[optind  + 2];
    }
  }

  return true;
}

bool HaveOption(const std::string& opt){
  return fldiagnostics_opts.count(opt);
}

const char* GetOption(const std::string& opt){
  assert(HaveOption(opt));

  int stat = 0;
  return GetOption(opt,  stat);
  }

const char* GetOption(const std::string& opt, int& stat){
  if(!HaveOption(opt)){
    stat = 1;
    return "";
  }

  stat = 0;
  return fldiagnostics_opts[opt].c_str();
}

void Error(const char* message){
  std::cerr << message << std::endl;
}

void FatalError(const char* message){
  Error(message);

  exit(1);
}

int main(int argc, char** argv){
#ifdef HAVE_MPI
  MPI::Init(argc, argv);
  // Undo some MPI init shenanigans
  chdir(getenv("PWD"));
#endif

  if(argc == 1){
    Help();
    return 0;
  }

  if(!ParseArgs(argc, argv)){
    Help();
    return 1;
  }

  std::queue< std::string > input_names, output_names;

  const char* input_name = GetOption("input_name");
  const char* output_name = GetOption("output_name");
  if(HaveOption("first_id")){
    for(unsigned int i = atoi(GetOption("first_id"));
      i <= (unsigned int)atoi(GetOption("last_id"));i++){
      char* buffer = (char*)malloc((i / 10 + 1) * sizeof(char));
      sprintf(buffer, "%u", i);
      input_names.push(std::string(input_name) + std::string("_")
        + std::string(buffer) + std::string(".vtu"));
      output_names.push(std::string(output_name)
        + std::string("_") + std::string(buffer) + std::string(".vtu"));
      free(buffer);
    }
  }
  else{
    input_names.push(std::string(input_name));
    output_names.push(std::string(output_name));
  }

  if(HaveOption("options_file")){
    std::string options_path;
    if(GetOption("options_file")[0] == '/'){
      options_path = std::string("");
    }
    else{
      options_path = std::string("/");
    }

    Spud::load_options((options_path + std::string(GetOption("options_file"))).c_str());
  }

  std::cout << "fldiagnostics\n"
            << "Command: " << GetOption("action") << "\n"
            << std::endl;

  while(input_names.size() > 0){
    const char* input_name = input_names.front().c_str();
      int input_name_len = strlen(input_name);
    const char* output_name = output_names.front().c_str();
      int output_name_len = strlen(output_name);

    const char* action = GetOption("action");
    if((strlen(action) == 8 and strncmp(action, "add-diag", 8) == 0) or (strlen(action) == 3 and strncmp(action, "add", 3) == 0)){
      const char* outfield_name = GetOption("outfield_name");
        int outfield_name_len = strlen(outfield_name);
      const char* meshfield_name = GetOption("meshfield_name");
        int meshfield_name_len = strlen(meshfield_name);
      const char* state_name;
        int state_name_len;

      if(HaveOption("state_name")){
        state_name = GetOption("state_name");
      }
      else{
        state_name = "";
      }
      state_name_len = strlen(state_name);
      int outfield_rank;
      if(HaveOption("outfield_rank")){
        outfield_rank = atoi(GetOption("outfield_rank"));
      }

      std::cout << "Adding diagnostic field: " << input_name << " + " << outfield_name << " => " << output_name << std::endl;

      fldiag_add_diag_fc(input_name, &input_name_len,
                         output_name, &output_name_len,
                         outfield_name, &outfield_name_len,
                         meshfield_name, &meshfield_name_len,
                         state_name, &state_name_len,
                         HaveOption("outfield_rank") ? &outfield_rank : NULL);
    }
/*
    else if(strlen(argv[1]) == 9 and strncmp(argv[1], "add-presc", 9) == 0){
      const char* outfield_name = GetOption("outfield_name");
        int outfield_name_len = strlen(outfield_name);
      const char* meshfield_name = GetOption("meshfield_name");
        int meshfield_name_len = strlen(meshfield_name);
      const char* state_name;
        int state_name_len;

      if(HaveOption("state_name")){
        state_name = GetOption("state_name");
      }
      else{
        state_name = "";
      }
      state_name_len = strlen(state_name);
      
      std::cout << "Adding prescribed field: " << input_name << " + " << outfield << " => " << output_name << std::endl;
      
      fldiag_add_presc_fc(input_name, &input_name_len,
                         output_name, &output_name_len,
                         outfield_name, &outfield_name_len,
                         meshfield_name, &meshfield_name_len,
                         state_name, &state_name_len)
    }
    else if(strlen(action) == 4 and strncmp(action, "diff", 4) == 0){
      const char* infield_name = GetOption("infield_name");
        int infield_name_len = strlen(infield_name);
      const char* outfield_name = GetOption("outfield_name");
        int outfield_name_len = strlen(outfield_name);
      const char* diff_type = GetOption("diff_type");
        int diff_type_len = strlen(diff_type);
      const char* coords_name = GetOption("coords_name");
        int coords_name_len = strlen(coords_name);

      fldiag_diff_fc(input_name, &input_name_len,
                     output_name, &output_name_len,
                     infield_name, &infield_name_len,
                     outfield_name, &outfield_name_len,
                     diff_type, &diff_type_len,
                     coords_name, &coords_name_len);
    }
*/
    else{
      std::cerr << "Command not found" << std::endl;
      return 1;
    }

    input_names.pop();
    output_names.pop();
  }

#ifdef HAVE_MPI
  MPI::Finalize();
#endif
  return 0;
}
