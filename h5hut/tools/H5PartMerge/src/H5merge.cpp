/* Description
 *
 *  This application merges and slices one, two or more H5hut files
 *  into a single one. The range of each file which should be merged/
 *  sliced can be specified on the command line where each input file
 *  is handled similar to a python array. Negative indices are counted
 *  from the end of the file i.e. -1 signifies the last step.
 *
 *  Copyright by Christof Kraus, 2007-2008, all rights reserved.
 *
 *  \author Christof Kraus
 *
 *  \date 2007 dec 27
 *
 *  \warning none
 *
 *  \attention none required
 *
 *  \bug this code is now ported to H5hut August 2011 Andreas Adelmann
 *
 *  \todo
 */

#include <hdf5.h>
#include "H5hut.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <cstdlib>

#include "optparse.hh"

#ifdef USE_BOOST
#include <boost/any.hpp>
#endif

#define MAX_LEN 128

using namespace std;

struct StringCmp {
    bool operator()(const string s1, const string s2) const {
        return s1 < s2;
    }
};

#ifdef USE_BOOST
using boost::any_cast;
typedef map<string, boost::any, StringCmp> StringAnyValueMap;
#endif

typedef map<string, int, StringCmp> StringIntMap;
typedef map<string, h5_int64_t, StringCmp> StringInt64Map;

h5_int64_t copyFileAttribute(h5_file_t *In, h5_file_t *Out, string AttrName,  h5_int64_t Type, h5_size_t Num);
h5_int64_t copyFileToStepAttribute(h5_file_t *In, h5_file_t *Out, string AttrName,  h5_int64_t Type, h5_size_t Num);

#ifdef USE_BOOST
h5_int64_t copyFileAttribute(h5_file_t *In, h5_file_t *Out, StringAnyValueMap &FileAttributes, string AttrName,
                             h5_int64_t Type, h5_size_t Num);
h5_int64_t copyFileToStepAttribute(h5_file_t *In, h5_file_t *Out, StringAnyValueMap &FileAttributes, string AttrName,
                                   h5_int64_t Type, h5_size_t Num);
#endif

h5_int64_t copyStepAttribute(h5_file_t *In, h5_file_t *Out, string AttrName,  h5_int64_t Type, h5_size_t Num);
h5_int64_t copyDataset(h5_file_t *In, h5_file_t *Out, string SetName,
                       h5_int64_t Type,  h5_float64_t *FloatArray, h5_int64_t *IntArray);

int main(int argc, char **argv) {

    using optparse::STORE;
    using optparse::STORE_TRUE;
    using optparse::STORE_FALSE;
    using optparse::STRING;
    using optparse::BOOL;
    using optparse::INT;
    using optparse::DOUBLE;

    h5_file_t *H5OutFile, *H5InFile;
    ifstream testopen;

    int NumInput = 0;
    int effNumInput = 0;
    string *InputFilenames = NULL;
    string Input;
    string Range;
    string from_str;
    int *from = NULL;
    string to_str;
    int *to = NULL;


    string OutputFilename;
    bool isOutputInput = false;
    int OutputIsInput = 0;

    char AttrName[MAX_LEN];
    int NumFileAttr;
    StringIntMap FileAttribName;
    StringInt64Map FileAttribType;
    StringInt64Map FileAttribNumType;

#ifdef USE_BOOST
    StringAnyValueMap FileAttributes;
#endif

    char SetName[MAX_LEN];
    int NumDataSets;
    StringIntMap DataSetName;
    StringInt64Map DataSetType;

    char StepAttrName[MAX_LEN];
    int NumStepAttr;
    StringIntMap StepAttribName;
    StringInt64Map StepAttribType;
    StringInt64Map StepAttribNumType;

    StringIntMap::iterator itName;
    StringInt64Map::iterator itType;
    StringInt64Map::iterator itNumType;


    h5_int64_t NumParticles = 0;
    h5_int64_t maxNumParticles = 0;

    h5_int64_t NumSteps;
    h5_int64_t effStep = 0;
    bool firstStepInput = true;

    h5_float64_t *FloatValues = NULL;
    h5_int64_t *IntegerValues = NULL;

    h5_int64_t Type;
    h5_int64_t SetType;
    h5_size_t NumType;

    h5_size_t SetNum;

    int verbosity = 0;
    h5_int64_t rc;

    string usage("H5merge\n\nmerges two or more H5hut files, but can also be used to slice a H5hut file\n\nUsage: H5merge [[options]] input1.h5[from1:to1] [input2.h5[from2:to2] ... inputN.h5[fromN:toN]] output.h5,\n       where fromX and toX are the step numbers which determine the range to be merged.\n       Negative values are counted from the end of the file\n\noptions:");
    optparse::OptionParser parser(usage);

    parser.add_option("-h", "--help", "help", "Show this message", STORE_TRUE, BOOL, "0");
    parser.add_option("-v", "--verbose", "verbosity", "increase verbosity", STORE, INT, "0");
    parser.add_option("-r", "--replace", "replace", "If output file name is the same as input file replace the input file at the end", STORE_TRUE, BOOL, "0");
    parser.add_option("-f", "--intersect-file-attributes", "intersect-file-attributes", "merge only the intersection of file attributes from all input files", STORE_TRUE, BOOL, "0");
    parser.add_option("-s", "--intersect-step-attributes", "intersect-step-attributes", "merge only the intersection of step attributes from all steps in all input files", STORE_TRUE, BOOL, "0");
    parser.add_option("-d", "--intersect-datasets", "intersect-datasets", "merge only the intersection of datasets of all input files", STORE_TRUE, BOOL, "0");
    try {
        parser.parse_args(argc, argv);
        if(parser.get<bool>("help")) {
            parser.help(cerr);
            exit(0);
        }
    } catch(optparse::OptionError e) {
        cerr << e.what() << endl << endl;
        parser.help(cerr);
        exit(1);
    }

    if(parser.arguments.size() < 2) {
        cerr << "No output file found." << endl << endl;
        parser.help(cerr);
        exit(1);
    }

    NumInput = parser.arguments.size(); //argc - 2;
    InputFilenames = new string[NumInput];
    from = new int[NumInput];
    to = new int[NumInput];

    verbosity = parser.get<int>("verbosity");
    H5SetVerbosityLevel(verbosity);

    /**************************************************************\
     * Parse input file names and ranges to be merged
    \**************************************************************/

    for(int i = 0; i < NumInput - 1; ++i) {
        Input = string(parser.arguments[i]);

        if(Input.find("[") != string::npos && Input.find("]") != string::npos) {
            Range = Input;
            Range.erase(0, Range.find_first_of("[") + 1);
            Range.erase(Range.find_first_of("]"));
            if(Range.find(":") != string::npos) {
                from_str = Range;
                from_str.erase(from_str.find_first_of(":"));
                if(from_str.length() == 0)
                    from[effNumInput] = 0;
                else
                    from[effNumInput] = atoi(from_str.c_str());

                to_str = Range;
                to_str.erase(0, to_str.find_first_of(":") + 1);
                if(to_str.length() == 0)
                    to[effNumInput] = -1;
                else
                    to[effNumInput] = atoi(to_str.c_str());
            } else {
                if(Range.length() > 0) {
                    from[effNumInput] = atoi(Range.c_str());
                    to[effNumInput] = from[effNumInput];
                } else {
                    from[effNumInput] = -1;
                    to[effNumInput] = 0;
                }
            }
            Input.erase(Input.find_first_of("["));

        } else {
            if(Input.find("[") != string::npos) {
                Input.erase(Input.find_first_of("["));
                cerr << "error when assigning a range for file " << Input << endl;
                exit(1);
            }
            if(Input.find("]") != string::npos) {
                Input.erase(Input.find_first_of("]"));
                cerr << "error when assigning a range for file " << Input << endl;
                exit(1);
            }
            from[effNumInput] = 0;
            to[effNumInput] = -1;
        }
        testopen.open(Input.c_str(), ios::in);
        if(testopen.good()) {
            InputFilenames[effNumInput] = Input;
            effNumInput++;
        } else {
            cerr << "Can't open " << Input << endl;
            exit(1);
        }
        testopen.close();
    }

    OutputFilename = string(parser.arguments[NumInput - 1]);

    /**************************************************************\
     * check whether output file name is also an input filename;
     * if true then append '_tmp' to the base name;
     * rename it in the end
    \**************************************************************/
    for(int i = 0; i < effNumInput; ++i) {
        if(InputFilenames[i] == OutputFilename) {
            if(!parser.get<bool>("replace")) {
                cout << "Do you really want to overwrite the input file " << InputFilenames[i] << "? [y/N] ";
                char answer = '0';
                cin.get(answer);

                if(answer != 'y' && answer != 'Y') {
                    cout << "Aborting..." << endl;
                    exit(1);
                }
            }
            ostringstream tmp;
            isOutputInput = true;
            OutputIsInput = i;
            OutputFilename.erase(OutputFilename.find_last_of("."));
            tmp << OutputFilename << "_tmp.h5";
            OutputFilename = tmp.str();
            break;
        }
    }

    /**************************************************************\
     * check whether output file allready exists
    \**************************************************************/
    if(!isOutputInput) {
        testopen.open(OutputFilename.c_str(), ios::in);
        if(testopen.good()) {
            testopen.close();
            if(!parser.get<bool>("replace")) {
                cout << "Do you really want to overwrite the file " << OutputFilename << "? [y/N] ";
                char answer = '0';
                cin.get(answer);

                if(answer != 'y' && answer != 'Y') {
                    cout << "Aborting..." << endl;
                    exit(1);
                }
            }
        } else if(parser.get<bool>("replace")) {
            cerr << "WARNING: option 'replace' chosen but output filename is not the same" << endl
                 << "         as any of the input filenames" << endl;
        }
    }

    for(int i = 0; i < effNumInput; ++i) {
        H5InFile = H5OpenFile(InputFilenames[i].c_str(), H5_O_RDONLY, 0);

        /**************************************************************\
         * fix range if values are negative
        \**************************************************************/
        int n_steps = H5GetNumSteps(H5InFile);

        if(from[i] >= n_steps)
            from[i] = n_steps - 1;
        else
            while(from[i] < 0)
                from[i] += n_steps;

        if(to[i] >= n_steps)
            to[i] = n_steps - 1;
        else
            while(to[i] < 0)
                to[i] += n_steps;
        if(to[i] < from[i]) {
            int tmp;
            tmp = to[i];
            to[i] = from[i];
            from[i] = tmp;
        }

        NumSteps = to[i] - from[i] + 1;

        for(int k = from[i]; k <= to[i]; ++k) {
            H5SetStep(H5InFile, k);

            /**************************************************************\
             * determine maximum number of particles in all input files
             * and all steps
            \**************************************************************/
            h5_int64_t localNumParticles = H5PartGetNumParticles(H5InFile);
            if(localNumParticles > maxNumParticles)
                maxNumParticles = localNumParticles;
        }


        /**************************************************************\
         * get the *intersection* of names of all file attributes
         * --------------- " -------------------  step attributes
         * --------------- " -------------------  datasets
         * from all files which are to be merged
        \**************************************************************/
        if(parser.get<bool>("intersect-file-attributes") ||
           parser.get<bool>("intersect-datasets")        ||
           parser.get<bool>("intersect-step-attributes")) {
            StringIntMap localStepAttribName;
            StringInt64Map localStepAttribType;
            StringInt64Map localStepAttribNumType;
            StringIntMap localDataSetName;
            StringInt64Map localDataSetType;
            StringIntMap::iterator itlocalName;
            StringInt64Map::iterator itlocalType;
            StringInt64Map::iterator itlocalNumType;


            if(parser.get<bool>("intersect-file-attributes")) {
                NumFileAttr = H5GetNumFileAttribs(H5InFile);
                for(int l = 0; l < NumFileAttr; ++l) {

                    H5GetFileAttribInfo(H5InFile, l, AttrName, MAX_LEN, &Type, &NumType);
                    itName = FileAttribName.find(string(AttrName));
                    if(itName == FileAttribName.end()) {
                        FileAttribName.insert(make_pair(string(AttrName), 1));
                        FileAttribType.insert(make_pair(string(AttrName), Type));
                        FileAttribNumType.insert(make_pair(string(AttrName), NumType));
                    } else {
                        if(FileAttribType[string(AttrName)] == Type) {
                            if(FileAttribNumType[string(AttrName)] == NumType) {
                                FileAttribName[string(AttrName)] += 1;
                            } else {
                                cerr << "WARNING: found same file attribute (" << AttrName << ") " << endl;
                                cerr << "         with same type but different length! Will be skipped!" << endl;
                            }
                        } else {
                            cerr << "WARNING: found same file attribute (" << AttrName << ") " << endl;
                            cerr << "         with different type! Will be skipped!" << endl;
                        }
                    }
                }
            }

            for(int k = from[i]; k <= to[i]; ++k) {
                H5SetStep(H5InFile, k);

                if(parser.get<bool>("intersect-datasets")) {
                    NumDataSets = H5PartGetNumDatasets(H5InFile);
                    for(h5_int64_t m = 0; m < NumDataSets; ++m) {
                        H5PartGetDatasetInfo(H5InFile, m, SetName, MAX_LEN, &SetType, &SetNum);
                        itName = localDataSetName.find(string(SetName));
                        if(itName == localDataSetName.end()) {
                            localDataSetName.insert(make_pair(string(SetName), 1));
                            localDataSetType.insert(make_pair(string(SetName), SetType));
                        } else {
                            if(localDataSetType[string(SetName)] == SetType) {
                                localDataSetName[string(SetName)] += 1;
                            } else {
                                cerr << "WARNING: found same dataset (" << SetName << ")" << endl;
                                cerr << "         with different type! Will be skipped!" << endl;
                            }
                        }
                    }
                }

                if(parser.get<bool>("intersect-step-attributes")) {
                    NumStepAttr = H5GetNumStepAttribs(H5InFile);
                    for(int m = 0; m < NumStepAttr; ++m) {
                        H5GetStepAttribInfo(H5InFile, m, StepAttrName, MAX_LEN, &Type, &NumType);
                        itlocalName = localStepAttribName.find(string(StepAttrName));
                        if(itlocalName == localStepAttribName.end()) {
                            localStepAttribName.insert(make_pair(string(StepAttrName), 1));
                            localStepAttribType.insert(make_pair(string(StepAttrName), Type));
                            localStepAttribNumType.insert(make_pair(string(StepAttrName), NumType));
                        } else {
                            if(localStepAttribType[string(StepAttrName)] == Type) {
                                if(localStepAttribNumType[string(StepAttrName)] == NumType) {
                                    localStepAttribName[string(StepAttrName)] += 1;
                                } else {
                                    cerr << "WARNING: found same step attribute (" << StepAttrName << ") " << endl;
                                    cerr << "         with same type but different length! Will be skipped!" << endl;
                                }
                            } else {
                                cerr << "WARNING: found same step attribute (" << StepAttrName << ")" << endl;
                                cerr << "         with different type! Will be skipped!" << endl;
                            }
                        }
                    }
                }
            }

            if(parser.get<bool>("intersect-datasets")) {
                // exclude those datasets which do not occur in all steps in file i
                for(itlocalName  = localDataSetName.begin(), itlocalType = localDataSetType.begin();
                    itlocalName != localDataSetName.end(); itlocalName++, itlocalType++) {
                    if(itlocalName->second != NumSteps) {
                        cerr << "WARNING: dataset '" << itlocalName->first << "' in file "
                             << InputFilenames[i] << " only found in " << itlocalName->second
                             << " of " << NumSteps << " steps;" << endl;
                        cerr << "         Skipping this dataset;" << endl;
                        localDataSetName.erase(itlocalName);
                    } else {
                        itName = DataSetName.find(itlocalName->first);
                        if(itName == DataSetName.end()) {
                            DataSetName.insert(make_pair(itlocalName->first, 1));
                            DataSetType.insert(make_pair(itlocalName->first, itlocalType->second));
                        } else {
                            if(DataSetType[itlocalName->first] == itlocalType->second) {
                                DataSetName[itlocalName->first] += 1;
                            } else {
                                cerr << "WARNING: found same dataset (" << itlocalName->first << ")" << endl;
                                cerr << "         with different type! Will be skipped!" << endl;
                            }
                        }
                    }
                }
            }

            if(parser.get<bool>("intersect-step-attributes")) {
                // exclude those step attributes which do not occur in all steps in file i
                itlocalType = localStepAttribType.begin();
                itlocalNumType = localStepAttribNumType.begin();
                for(itlocalName = localStepAttribName.begin(); itlocalName != localStepAttribName.end();
                    itlocalName++, itlocalType++, itlocalNumType++) {
                    if(itlocalName->second != NumSteps) {
                        cerr << "WARNING: step attribute '" << itlocalName->first << "' in file " << InputFilenames[i]
                             << " only found in " << itlocalName->second << " of " << NumSteps << " steps;" << endl;
                        cerr << "         Skipping this step attribute;" << endl;
                    } else {
                        itName = StepAttribName.find(itlocalName->first);
                        if(itName == StepAttribName.end()) {
                            StepAttribName.insert(make_pair(itlocalName->first, 1));
                            StepAttribType.insert(make_pair(itlocalName->first, itlocalType->second));
                            StepAttribNumType.insert(make_pair(itlocalName->first, itlocalNumType->second));
                        } else {
                            if(StepAttribType[itlocalName->first] == itlocalType->second) {
                                if(StepAttribNumType[itlocalName->first] == itlocalNumType->second) {
                                    StepAttribName[itlocalName->first] += 1;
                                } else {
                                    cerr << "WARNING: found same step attribute (" << itlocalName->first << ") " << endl;
                                    cerr << "         with same type but different length! Will be skipped!" << endl;
                                }
                            } else {
                                cerr << "WARNING: found same step attribute (" << itlocalName->first << ")" << endl;
                                cerr << "         with different type! Will be skipped!" << endl;
                            }
                        }
                    }
                }
            }

            H5CloseFile(H5InFile);
        }
    }


    if(parser.get<bool>("intersect-file-attributes")) {
        //erase all those file attributes which do not occur in all files
        itType = FileAttribType.begin();
        itNumType = FileAttribNumType.begin();
        for(itName = FileAttribName.begin(); itName != FileAttribName.end(); itName++, itType++, itNumType++) {
            if(itName->second != effNumInput) {
                cerr << "WARNING: file attribute '" << itName->first << "' only found in " << itName->second << " of "
                     << effNumInput << " input files;" << endl;
                cerr << "         skipping this attribute!" << endl;
                FileAttribName.erase(itName);
                FileAttribType.erase(itType);
                FileAttribNumType.erase(itNumType);
            }
        }

        if(verbosity > 0) {
            cout << "found the following file attributes (intersection of all file attributes):" << endl;
            for(itName = FileAttribName.begin(); itName != FileAttribName.end(); ++itName) {
                cout << "file attribute '" << itName->first << "'" << endl;
            }
            cout << endl;
        }
    }

    if(parser.get<bool>("intersect-datasets")) {
        //erase all those datasets which do not occur in all steps of all files
        for(itName = DataSetName.begin(), itType = DataSetType.begin(); itName != DataSetName.end(); itName++, itType++) {
            if(itName->second != effNumInput) {
                cerr << "WARNING: dataset '" << itName->first << "' only found in " << itName->second
                     << " of " << effNumInput << " input files;" << endl;
                cerr << "         skipping this data set!" << endl;
                DataSetName.erase(itName);
                DataSetType.erase(itType);
            }
        }

        if(verbosity > 0) {
            cout << "found the following data sets (intersection of all datasets):" << endl;
            for(itName = DataSetName.begin(); itName != DataSetName.end(); itName++)
                cout << "dataset '" << itName->first << "'" << endl;
            cout << endl;
        }
    }

    if(parser.get<bool>("intersect-step-attributes")) {
        //erase all those step attributes which do not occur in all steps of all fiels
        itType = StepAttribType.begin();
        itNumType = StepAttribNumType.begin();
        for(itName = StepAttribName.begin(); itName != StepAttribName.end(); itName++, itType++, itNumType++) {
            if(itName->second != effNumInput) {
                cerr << "WARNING: step attribute '" << itName->first << "' only found in " << itName->second
                     << " of " << effNumInput << " input files;" << endl;
                cerr << "         skipping this attribute!" << endl;
                StepAttribName.erase(itName);
                StepAttribType.erase(itType);
                StepAttribNumType.erase(itNumType);
            }
        }

        if(verbosity > 0) {
            cout << "found the following step attributes (intersection of all step attributes):" << endl;
            for(itName = StepAttribName.begin(); itName != StepAttribName.end(); ++itName) {
                cout << "step attribute '" << itName->first << "'" << endl;
            }
            cout << endl;
        }
    }


    H5OutFile = H5OpenFile(OutputFilename.c_str(), H5_O_WRONLY, 0);
    if(H5OutFile == NULL) {
        cerr << "ABORT: could not open '" << OutputFilename << "' for writing!" << endl;
        exit(1);
    }

    FloatValues = new h5_float64_t[maxNumParticles];
    IntegerValues = new h5_int64_t[maxNumParticles];

    for(int i = 0; i < effNumInput; ++i) {
        H5InFile = H5OpenFile(InputFilenames[i].c_str(), H5_O_RDONLY, 0);

        H5SetStep(H5InFile, from[i]);
        H5SetStep(H5OutFile, effStep);

        if(parser.get<bool>("intersect-file-attributes")) {
            itType = FileAttribType.begin();
            itNumType = FileAttribNumType.begin();
            for(itName = FileAttribName.begin(); itName != FileAttribName.end(); itName++, itType++, itNumType++) {
                if(i == 0) {
#ifdef USE_BOOST
                    rc = copyFileAttribute(H5InFile, H5OutFile, FileAttributes, itName->first, itType->second, itNumType->second);
#else
                    rc = copyFileAttribute(H5InFile, H5OutFile, itName->first, itType->second, itNumType->second);
#endif
                } else {
#ifdef USE_BOOST
                    rc = copyFileToStepAttribute(H5InFile, H5OutFile, FileAttributes, itName->first, itType->second, itNumType->second);
#else
                    rc = copyFileToStepAttribute(H5InFile, H5OutFile, itName->first, itType->second, itNumType->second);
#endif
                }
                if(rc != H5_SUCCESS)
                    cerr << "WARNING: could not write file attribute '" << AttrName << "'" << endl;
            }
        } else {
            h5_int64_t NumFileAttr = H5GetNumFileAttribs(H5InFile);
            for(h5_int64_t l = 0; l < NumFileAttr; ++l) {
                H5GetFileAttribInfo(H5InFile, l, AttrName, MAX_LEN, &Type, &NumType);
                if(i == 0) {
#ifdef USE_BOOST
                    rc = copyFileAttribute(H5InFile, H5OutFile, FileAttributes, string(AttrName), Type, NumType);
#else
                    rc = copyFileAttribute(H5InFile, H5OutFile, string(AttrName), Type, NumType);
#endif
                } else {
#ifdef USE_BOOST
                    rc = copyFileToStepAttribute(H5InFile, H5OutFile, FileAttributes, string(AttrName), Type, NumType);
#else
                    rc = copyFileToStepAttribute(H5InFile, H5OutFile, string(AttrName), Type, NumType);
#endif
                }
                if(rc != H5_SUCCESS)
                    cerr << "WARNING: could not write file attribute '" << AttrName << "'" << endl;
            }
        }

        firstStepInput = true;
        for(int l = from[i]; l <= to[i]; ++l, ++effStep) {
            if(!firstStepInput) {
                rc = H5SetStep(H5InFile, l);
                if(rc != H5_SUCCESS) {
                    cerr << "ABORT: could not change to step #" << l << " in input file '"
                         << InputFilenames[i] << "'" << endl;
                    rc = H5CloseFile(H5InFile);
                    rc = H5CloseFile(H5OutFile);
                    exit(1);
                }
                rc = H5SetStep(H5OutFile, effStep);
                if(rc != H5_SUCCESS) {
                    cerr << "ABORT: could not change to step #" << effStep << " in output file '"
                         << OutputFilename << "'" << endl;
                    rc = H5CloseFile(H5InFile);
                    rc = H5CloseFile(H5OutFile);
                    exit(1);
                }
            } else {
                firstStepInput = false;
            }
            NumParticles = H5PartGetNumParticles(H5InFile);
            rc = H5PartSetNumParticles(H5OutFile, NumParticles);
            if(rc != H5_SUCCESS) {
                cerr << "ABORT: could not set the number of particles in output file '" << OutputFilename << "'" << endl;
                rc = H5CloseFile(H5InFile);
                rc = H5CloseFile(H5OutFile);
                exit(1);
            }
            if(parser.get<bool>("intersect-step-attributes")) {
                itType = StepAttribType.begin();
                itNumType = StepAttribNumType.begin();
                for(itName = StepAttribName.begin(); itName != StepAttribName.end(); ++itName, ++itType, ++itNumType) {
                    rc = copyStepAttribute(H5InFile, H5OutFile, itName->first, itType->second, itNumType->second);
                    if(rc != H5_SUCCESS)
                        cerr << "WARNING: could not write step attribute '" << itName->first << "'" << endl;
                }
            } else {
                h5_int64_t NumStepAttr = H5GetNumStepAttribs(H5InFile);
                for(h5_int64_t m = 0; m < NumStepAttr; ++m) {
                    H5GetStepAttribInfo(H5InFile, m, StepAttrName, MAX_LEN, &Type, &NumType);
                    rc = copyStepAttribute(H5InFile, H5OutFile, string(StepAttrName), Type, NumType);
                    if(rc != H5_SUCCESS)
                        cerr << "WARNING: could not write step attribute '" << StepAttrName << "'" << endl;
                }
            }

            if(parser.get<bool>("intersect-datasets")) {

                for(itType = DataSetType.begin(); itType != DataSetType.end(); ++itType) {
                    rc = copyDataset(H5InFile, H5OutFile, itType->first, itType->second, FloatValues, IntegerValues);
                    if(rc != H5_SUCCESS)
                        cerr << "WARNING: could not write dataset '" << itType->first << "'" << endl;
                }
            } else {
                h5_int64_t NumDataSets = H5PartGetNumDatasets(H5InFile);
                for(h5_int64_t m = 0; m < NumDataSets; ++m) {
                    H5PartGetDatasetInfo(H5InFile, m, SetName, MAX_LEN, &SetType, &SetNum);
                    rc = copyDataset(H5InFile, H5OutFile, string(SetName), SetType, FloatValues, IntegerValues);
                    if(rc != H5_SUCCESS)
                        cerr << "WARNING: could not write dataset '" << SetName << "'" << endl;
                }
            }

        }
        rc = H5CloseFile(H5InFile);
    }
    rc = H5CloseFile(H5OutFile);

    if(isOutputInput) {
        rename(OutputFilename.c_str(), InputFilenames[OutputIsInput].c_str());
    }

    delete[] FloatValues;
    delete[] IntegerValues;
    delete[] InputFilenames;
    delete[] from;
    delete[] to;
}

h5_int64_t copyFileAttribute(h5_file_t *In, h5_file_t *Out, string AttrName,  h5_int64_t Type, h5_size_t Num) {
    h5_int64_t rc;
    if(Type == H5T_NATIVE_DOUBLE) {
        void *Attrib = (double *)malloc(sizeof(double) * Num);
        rc = H5ReadFileAttribFloat64(In, AttrName.c_str(), (double *)Attrib);
        if(rc == H5_SUCCESS)
            rc = H5WriteFileAttribFloat64(Out, AttrName.c_str(), (double *)Attrib, Num);
        free(Attrib);
        return rc;
    } else if(Type == H5T_NATIVE_CHAR) {
        void *Attrib = (char *)malloc(sizeof(char) * Num);
        rc = H5ReadFileAttribString(In, AttrName.c_str(), (char *)Attrib);
        if(rc == H5_SUCCESS)
            rc = H5WriteFileAttribString(Out, AttrName.c_str(), (char *)Attrib);
        free(Attrib);
        return rc;
    } else if(Type == H5T_NATIVE_INT64) {
        void *Attrib = (h5_int64_t *)malloc(sizeof(h5_int64_t) * Num);
        rc = H5ReadFileAttribInt64(In, AttrName.c_str(), (h5_int64_t *)Attrib);
        if(rc == H5_SUCCESS)
            rc = H5WriteFileAttribInt64(Out, AttrName.c_str(), (h5_int64_t *)Attrib, Num);
        free(Attrib);
        return rc;
    } else {
        cerr << "WARNING: unknown type for file attribue " << AttrName << "! Skipping it. " << endl;
        return H5_ERR_INVAL;
    }
}


h5_int64_t copyFileToStepAttribute(h5_file_t *In, h5_file_t *Out, string AttrName,  h5_int64_t Type, h5_size_t Num) {
    h5_int64_t rc;
    if(Type == H5T_NATIVE_DOUBLE) {
        void *Attrib = (double *)malloc(sizeof(double) * Num);
        rc = H5ReadFileAttribFloat64(In, AttrName.c_str(), (double *)Attrib);
        if(rc == H5_SUCCESS)
            rc = H5WriteStepAttribFloat64(Out, AttrName.c_str(), (double *)Attrib, Num);
        free(Attrib);
        return rc;
    } else if(Type == H5T_NATIVE_CHAR) {
        void *Attrib = (char *)malloc(sizeof(char) * Num);
        rc = H5ReadFileAttribString(In, AttrName.c_str(), (char *)Attrib);
        if(rc == H5_SUCCESS)
            rc = H5WriteStepAttribString(Out, AttrName.c_str(), (char *)Attrib);
        free(Attrib);
        return rc;
    } else if(Type == H5T_NATIVE_INT64) {
        void *Attrib = (h5_int64_t *)malloc(sizeof(h5_int64_t) * Num);
        rc = H5ReadFileAttribInt64(In, AttrName.c_str(), (h5_int64_t *)Attrib);
        if(rc == H5_SUCCESS)
            rc = H5WriteStepAttribInt64(Out, AttrName.c_str(), (h5_int64_t *)Attrib, Num);
        free(Attrib);
        return rc;
    } else {
        cerr << "WARNING: unknown type for file attribue " << AttrName << "! Skipping it. " << endl;
        return H5_ERR_INVAL;
    }
}

h5_int64_t copyStepAttribute(h5_file_t *In, h5_file_t *Out, string AttrName,  h5_int64_t Type, h5_size_t Num) {
    h5_int64_t rc;
    if(Type == H5T_NATIVE_DOUBLE) {
        void *Attrib = (double *)malloc(sizeof(double) * Num);;
        rc = H5ReadStepAttribFloat64(In, AttrName.c_str(), (double *)Attrib);
        if(rc == H5_SUCCESS)
            rc = H5WriteStepAttribFloat64(Out, AttrName.c_str(), (double *)Attrib, Num);
        free(Attrib);
        return rc;
    } else if(Type == H5T_NATIVE_CHAR) {
        void *Attrib = (char *)malloc(sizeof(char) * Num);;
        rc = H5ReadStepAttribString(In, AttrName.c_str(), (char *)Attrib);
        if(rc == H5_SUCCESS)
            rc = H5WriteStepAttribString(Out, AttrName.c_str(), (char *)Attrib);
        free(Attrib);
        return rc;
    } else if(Type == H5T_NATIVE_INT64) {
        void *Attrib = (h5_int64_t *)malloc(sizeof(h5_int64_t) * Num);
        rc = H5ReadStepAttribInt64(In, AttrName.c_str(), (h5_int64_t *)Attrib);
        if(rc == H5_SUCCESS)
            rc = H5WriteStepAttribInt64(Out, AttrName.c_str(), (h5_int64_t *)Attrib, Num);
        free(Attrib);
        return rc;
    } else {
        cerr << "WARNING: unknown type for step attribute " << AttrName << "! Skipping it." << endl;
        return H5_ERR_INVAL;
    }
}


h5_int64_t copyDataset(h5_file_t *In, h5_file_t *Out, string SetName, h5_int64_t Type,  h5_float64_t *FloatArray, h5_int64_t *IntegerArray) {
    h5_int64_t rc;
    if(Type == H5T_NATIVE_DOUBLE) {
        rc = H5PartReadDataFloat64(In, SetName.c_str(), FloatArray);
        if(rc == H5_SUCCESS)
            rc = H5PartWriteDataFloat64(Out, SetName.c_str(), FloatArray);
        return rc;
    } else if(Type == H5T_NATIVE_INT64) {
        rc = H5PartReadDataInt64(In, SetName.c_str(), IntegerArray);
        if(rc == H5_SUCCESS)
            rc = H5PartWriteDataInt64(Out, SetName.c_str(), IntegerArray);
        return rc;
    } else {
        cerr << "WARNING: unknown type for dataset " << SetName << "! Skipping it." << endl;
        return H5_ERR_INVAL;
    }
}

#ifdef USE_BOOST
h5_int64_t copyFileAttribute(h5_file_t *In, h5_file_t *Out, StringAnyValueMap &FileAttributes, string AttrName,  h5_int64_t Type, h5_int64_t Num) {
    h5_int64_t rc;
    if(Type == H5T_NATIVE_DOUBLE) {
        void *Attrib = (double *)malloc(sizeof(double) * Num);
        rc = H5PartReadFileAttrib(In, AttrName.c_str(), (double *)Attrib);
        if(rc == H5_SUCCESS)
            rc = H5PartWriteFileAttrib(Out, AttrName.c_str(), Type, (double *)Attrib, Num);
        if(rc == H5_SUCCESS)
            FileAttributes.insert(make_pair(AttrName, boost::any(*(double *)Attrib)));
        free(Attrib);
        return rc;
    } else if(Type == H5T_NATIVE_CHAR) {
        void *Attrib = (char *)malloc(sizeof(char) * Num);
        rc = H5PartReadFileAttrib(In, AttrName.c_str(), (char *)Attrib);
        if(rc == H5_SUCCESS)
            rc = H5PartWriteFileAttrib(Out, AttrName.c_str(), Type, (char *)Attrib, Num);
        if(rc == H5_SUCCESS)
            FileAttributes.insert(make_pair(AttrName, boost::any(string((char *)Attrib))));
        free(Attrib);
        return rc;
    } else if(Type == H5T_NATIVE_INT64) {
        void *Attrib = (h5_int64_t *)malloc(sizeof(h5_int64_t) * Num);
        rc = H5PartReadFileAttrib(In, AttrName.c_str(), (h5_int64_t *)Attrib);
        if(rc == H5_SUCCESS)
            rc = H5PartWriteFileAttrib(Out, AttrName.c_str(), Type, (h5_int64_t *)Attrib, Num);
        if(rc == H5_SUCCESS)
            FileAttributes.insert(make_pair(AttrName, boost::any(*(h5_int64_t *)Attrib)));
        free(Attrib);
        return rc;
    } else {
        cerr << "WARNING: unknown type for file attribue " << AttrName << "! Skipping it. " << endl;
        return H5_ERR_INVAL;
    }
}


h5_int64_t copyFileToStepAttribute(h5_file_t *In, h5_file_t *Out, StringAnyValueMap &FileAttributes, string AttrName,  h5_int64_t Type, h5_int64_t Num) {
    h5_int64_t rc;
    if(Type == H5T_NATIVE_DOUBLE) {
        void *Attrib = (double *)malloc(sizeof(double) * Num);
        rc = H5PartReadFileAttrib(In, AttrName.c_str(), (double *)Attrib);
        if(rc == H5_SUCCESS) {
            if(FileAttributes.find(AttrName) == FileAttributes.end()) {
                rc = H5PartWriteStepAttrib(Out, AttrName.c_str(), Type, (double *)Attrib, Num);
                if(rc == H5_SUCCESS)
                    FileAttributes.insert(make_pair(AttrName, boost::any(*(double *)Attrib)));
            } else if(any_cast<double>(FileAttributes[AttrName]) != *(double *)Attrib) {
                rc = H5PartWriteStepAttrib(Out, (AttrName + " changed to").c_str(), Type, (double *)Attrib, Num);
                if(rc == H5_SUCCESS)
                    FileAttributes[AttrName] = boost::any(*(double *)Attrib);
            }
        }
        free(Attrib);
        return rc;
    } else if(Type == H5T_NATIVE_CHAR) {
        void *Attrib = (char *)malloc(sizeof(char) * Num);
        rc = H5PartReadFileAttrib(In, AttrName.c_str(), (char *)Attrib);
        if(rc == H5_SUCCESS) {
            if(FileAttributes.find(AttrName) == FileAttributes.end()) {
                rc = H5PartWriteStepAttrib(Out, AttrName.c_str(), Type, (char *)Attrib, Num);
                if(rc == H5_SUCCESS)
                    FileAttributes.insert(make_pair(AttrName, boost::any(string((char *)Attrib))));
            } else if(strcmp((any_cast<string>(FileAttributes[AttrName])).c_str(), (char *)Attrib) != 0) {
                rc = H5PartWriteStepAttrib(Out, (AttrName + " changed to").c_str(), Type, (char *)Attrib, Num);
                if(rc == H5_SUCCESS)
                    FileAttributes[AttrName] = boost::any(string((char *)Attrib));
            }
        }
        free(Attrib);
        return rc;
    } else if(Type == H5T_NATIVE_INT64) {
        void *Attrib = (h5_int64_t *)malloc(sizeof(h5_int64_t) * Num);
        rc = H5PartReadFileAttrib(In, AttrName.c_str(), (h5_int64_t *)Attrib);
        if(rc == H5_SUCCESS) {
            if(FileAttributes.find(AttrName) == FileAttributes.end()) {
                rc = H5PartWriteStepAttrib(Out, AttrName.c_str(), Type, (h5_int64_t *)Attrib, Num);
                FileAttributes.insert(make_pair(AttrName, boost::any((h5_int64_t *)Attrib)));
            } else if(any_cast<h5_int64_t>(FileAttributes[AttrName]) != *(h5_int64_t *)Attrib) {
                rc = H5PartWriteStepAttrib(Out, (AttrName + " changed to").c_str(), Type, (h5_int64_t *)Attrib, Num);
                if(rc == H5_SUCCESS)
                    FileAttributes[AttrName] = boost::any((h5_int64_t *)Attrib);
            }
        }
        free(Attrib);
        return rc;
    } else {
        cerr << "WARNING: unknown type for file attribue " << AttrName << "! Skipping it. " << endl;
        return H5_ERR_INVAL;
    }
}
#endif
