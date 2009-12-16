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

#include <confdefs.h>

#include <vtk.h>

#include <cmath>
#include <iostream>
#include <cassert>
#include <getopt.h>
#include <string>
#include <map>
#include <vector>

using namespace std;

// the function that actually does the work - in the femtools dir
extern int projections(int nPoints, double *x, double *y, double *z, string current_coord, string output_coord);


#define DEBUG(msg) if (OptionIsSet("verbose")) cerr << msg << endl;

void PrintUsageInformation(const char* executableName);

enum ExitCodes{
  SUCCESS = 0,
  BAD_ARGUMENTS = -1,
  UNKNOWN_CELL_TYPE = -2
};

map<string, string> optionMap;

bool OptionIsSet(const string& option){
  return optionMap.find(option) != optionMap.end();
}

void ParseArguments(int argc, char** argv){
  struct option longOptions[] = {
    {"help",       0, 0, 'h'},
    {"incoords",   1, 0, 'I'},
    {"outcoords",  1, 0, 'O'},
    {"output",     1, 0, 'o'},
    {"verbose",    0, 0, 'v'},
    {0,            0, 0, 0}
  };
  
  int optionIndex = 0;
  int c;
  
  // Make sure we have enough options
  if (argc>=5){    
    optionMap["infile"] = argv[argc-1];
    argc--;

    bool flag = false;
    fstream fin;
    fin.open(optionMap["infile"].c_str(),ios::in);
    if(!fin.is_open()){
      cerr<<"ERROR: no such file: "<<optionMap["infile"]<<endl;
      flag=true;
    }
    fin.close();
  }
  
  // Parse remaining options
  while (true){
    c = getopt_long(argc, argv, "h:I:O:o:v", longOptions, &optionIndex);
    if (c == -1) break;
    
    switch (c){  
    case 'h':
      PrintUsageInformation(argv[0]);
      exit(SUCCESS);

    case 'I':
      optionMap["incoord"] = optarg;
      break;

    case 'O':
      optionMap["outcoord"] = optarg;
      break;
      
    case 'o':
      optionMap["outfile"] = optarg;
      break;
            
    case 'v': 
      optionMap["verbose"] = "true";
      break;
            
    case '?':
      PrintUsageInformation(argv[0]);
      exit(BAD_ARGUMENTS);
    }
  }
  
  // Make sure we have the right options
  if (!OptionIsSet("incoord")){
    cout << "Coordinate type of input file not specified!" << endl;
    PrintUsageInformation(argv[0]);
    exit(BAD_ARGUMENTS);
  }
  if (!OptionIsSet("outcoord")){
    cout << "Coordinate type of output not specified!" << endl;
    PrintUsageInformation(argv[0]);
    exit(BAD_ARGUMENTS);
  }
}

void PrintUsageInformation(const char* executableName){
  DEBUG("void PrintUsageInformation(const char* executableName)");
  cerr << "Usage: " << executableName << "[OPTIONS] -I <in-coordinates> -O <out-coordinates> infile.vtu" << endl << endl
    
       << "Options:" << endl

       << "-h, --help" << endl
       << "\tPrint this usage information" << endl << endl

       << "-I, --incoords" << endl
       << "\tCoordinates of input file. Valid types are:" << endl
       << "\t type      | description "<< endl
       << "\t-------------------------"<< endl
       << "\t cart      | Cartesian (meters)" << endl
       << "\t spherical | Longitude/Latitude" << endl
       << "\t stereo    | Stereographic Projection" << endl << endl

       << "-O, --outcoords" << endl
       << "\tCoordinates of output file. Valid types are:" << endl
       << "\t type      | description "<< endl
       << "\t-------------------------"<< endl
       << "\t cart      | Cartesian (meters)" << endl
       << "\t spherical | Longitude/Latitude" << endl
       << "\t stereo    | Stereographic Projection" << endl << endl

       << "-o, --output=FILENAME.vtu" << endl
       << "\tFile that the result is outputed to. **The default is to overwrite the input file**" << endl << endl
    
       << "-v, --verbose" << endl
       << "\tBe verbose" << endl;
}

vtkUnstructuredGrid* ReadGrid(const string& filename){
  DEBUG("vtkUnstructuredGrid* ReadGrid("<<filename<<")");
  vtkXMLUnstructuredGridReader* reader = vtkXMLUnstructuredGridReader::New();
  
  reader->SetFileName(filename.c_str());
  reader->Update();
  
  vtkUnstructuredGrid* grid = vtkUnstructuredGrid::New();
  grid->DeepCopy(reader->GetOutput());
  grid->Update();

  reader->Delete();
  return grid;
}


void WriteGrid(vtkUnstructuredGrid* grid, const string& filename){
  DEBUG("void WriteGrid(vtkUnstructuredGrid* grid, const string& filename)");
  vtkXMLUnstructuredGridWriter* writer= vtkXMLUnstructuredGridWriter::New();
  writer->SetFileName( filename.c_str() );
  writer->SetInput(grid);
  vtkZLibDataCompressor* compressor = vtkZLibDataCompressor::New();
  writer->SetCompressor(compressor);
  writer->Write();
  writer->Delete();
  compressor->Delete();  
}

int main(int argc, char** argv){
  ParseArguments(argc, argv);
  
  // Load
  vtkUnstructuredGrid* grid = ReadGrid( optionMap["infile"] );
  
  // Point definitions
  vtkIdType npts = grid->GetNumberOfPoints();
  cerr<<"num of points - "<<npts<<endl;

  double *x = new double[npts];
  double *y = new double[npts];
  double *z = new double[npts];
  
  for(vtkIdType i=0;i<npts;i++){
    double xyz[3];
    grid->GetPoint(i, xyz);
    x[i] = xyz[0];
    y[i] = xyz[1];
    z[i] = xyz[2];
  }

  cout<<"Converting from : "<<optionMap["incoord"]<<" to: "<<optionMap["outcoord"]<<endl;

  int err = projections(npts, x, y, z, optionMap["incoord"], optionMap["outcoord"]);
   
  for(vtkIdType i=0;i<npts;i++){
    grid->GetPoints()->SetPoint(i, x[i], y[i], z[i]);
  }

  // Write
  if ( OptionIsSet("outfile") ){
    WriteGrid( grid, optionMap["outfile"] );
  }else{
    WriteGrid( grid, optionMap["infile"] );
  }

  return SUCCESS;
}
