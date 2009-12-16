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
#include <set>
#include <string>
#include <map>
#include <vector>

#define DEBUG(msg) if (OptionIsSet("debug")) cout << msg << endl;


using namespace std;

map<string, string> optionMap;

const double PI = 3.1415926535897932;
const double PIx4 = 4 * PI;

enum ExitCodes{
  SUCCESS = 0,
  BAD_ARGUMENTS,
  UNKNOWN_CELL_TYPE
};

struct Vector3d{
  double x, y, z;
  
  Vector3d(): x(0), y(0), z(0) {}
  
  Vector3d(const double xyz[3]): x(xyz[0]), y(xyz[1]), z(xyz[2]) {}
  
  double Magnitude() const{
    return sqrt(x*x + y*y + z*z);
  }
  
  void Zero(){
    x = y = z = 0;
  }
  
  Vector3d& Multiply(const Vector3d& v, double s){
    x = v.x * s;
    y = v.y * s;
    z = v.z * s;
    
    return *this;
  }
  
  Vector3d& Add(const Vector3d& v1, const Vector3d& v2){
    x = v1.x + v2.x;
    y = v1.y + v2.y;
    z = v1.z + v2.z;
    
    return *this;
  }
};

ostream& operator<<(ostream& stream, const Vector3d& v){
  stream << "(" << v.x << ", " << v.y << ", " << v.z << ")";
  
  return stream;
}

bool OptionIsSet(const string& option){
  return optionMap.find(option) != optionMap.end();
}

vtkUnstructuredGrid* LoadGrid(const string& filename){
  vtkXMLUnstructuredGridReader* reader = vtkXMLUnstructuredGridReader::New();
  
  reader->SetFileName( filename.c_str() );
  reader->Update();
  
  vtkUnstructuredGrid* grid = vtkUnstructuredGrid::New();
  grid->DeepCopy(reader->GetOutput());
  grid->Update();
  
  reader->Delete();
  
  return grid;
}

void WriteGrid(vtkUnstructuredGrid* grid, const string& filename){
  if (OptionIsSet("radial-scaling")){
    double rscale=strtod(optionMap["radial-scaling"].c_str(), NULL);
    
    vtkIdType npts = grid->GetNumberOfPoints();
    double xyz[3];
    const double earth_radius=6378000.0;
    
    vtkDoubleArray* depthv = vtkDoubleArray::New();
    depthv->SetNumberOfComponents(1);
    depthv->SetNumberOfTuples(npts);  
    depthv->SetName("depth");
    
    for(vtkIdType i=0;i<npts; i++){
      // Deal with verticies
      grid->GetPoints()->GetPoint(i, xyz);
      
      // Convert to spherical polar coordinates
      double r = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]);
      double longitude = atan2(xyz[1], xyz[0]);
      double colatitude = acos(xyz[2]/r);
      
      // Scale along r
      double depth = r - earth_radius;
      depthv->SetTuple1(i, depth);
      
      depth*=rscale;
      
      r = earth_radius+depth;
      
      // Convert back to Cartesian
      double new_xyz[3];
      new_xyz[0] = r*sin(colatitude)*cos(longitude);
      new_xyz[1] = r*sin(colatitude)*sin(longitude);
      new_xyz[2] = r*cos(colatitude);
      
      grid->GetPoints()->SetPoint(i, new_xyz);
      
      if(grid->GetPointData()->GetVectors()!=NULL){ // Deal with vectors
  double v[3];
  grid->GetPointData()->GetVectors()->GetTuple(i, v);
  
  for(size_t j=0;j<3;j++)
    v[j] += xyz[j];
  
  r=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  if(r>0.0){
    longitude = atan2(v[1], v[0]);
    colatitude = acos(v[2]/r);
    
    // Scale along r 
    depth = r - earth_radius;
    depth*=rscale;
    
    r = earth_radius+depth;
    
    // Convert back to Cartesian
    double new_v[3];
    new_v[0] = r*sin(colatitude)*cos(longitude);
    new_v[1] = r*sin(colatitude)*sin(longitude);
    new_v[2] = r*cos(colatitude);
    
    for(size_t j=0;j<3;j++)
      new_v[j] -= new_xyz[j];
    
    grid->GetPointData()->GetVectors()->SetTuple(i, new_v);
  }
      }
    }
    grid->GetPointData()->AddArray(depthv);
  }
  
  // -  
  vector< set<int> > graph(grid->GetNumberOfPoints());
  for(vtkIdType i=0;i<grid->GetNumberOfCells();i++){
    vtkCell *cell = grid->GetCell(i);
    
    for(int j=0; j<cell->GetNumberOfPoints(); j++){
      vtkIdType pt0 = cell->GetPointId(j);
      for(vtkIdType k=0; k<cell->GetNumberOfPoints(); k++){
  graph[pt0].insert(cell->GetPointId(k));
      }
    }
  }
  
  vector<int> distance(grid->GetNumberOfPoints(), -1);
  set<int> front, new_front;
  front.insert(0);
  int level=0;
  distance[0] = level;
 
  while(!front.empty()){
    level++;
    for(set<int>::const_iterator it=front.begin();it!=front.end(); it++){
      for(set<int>::const_iterator jt=graph[*it].begin();jt!=graph[*it].end();jt++){
  if(distance[*jt]==-1){
    new_front.insert(*jt);
    distance[*jt]=level;
  }
      }
    }

    front.swap(new_front);
    new_front.clear();
  }
  
  vtkDoubleArray* graph_distance = vtkDoubleArray::New();
  graph_distance->SetNumberOfComponents(1);
  graph_distance->SetNumberOfTuples(grid->GetNumberOfPoints());  
  graph_distance->SetName("Graph distance");
  
  for(int i=0;i<grid->GetNumberOfPoints();i++){
    double flt=distance[i];
    graph_distance->SetTuple(i, &flt);
  }

  grid->GetPointData()->AddArray(graph_distance);

  vtkXMLUnstructuredGridWriter* writer= vtkXMLUnstructuredGridWriter::New();
  writer->SetFileName( filename.c_str() );
  writer->SetInput(grid);
  vtkZLibDataCompressor* compressor = vtkZLibDataCompressor::New();
  writer->SetCompressor(compressor);
  writer->Write();
  writer->Delete();
  compressor->Delete();  
}

double CalculateCellVolume(vtkCell *cell){
  double r0[3], r1[3], r2[3], r3[3], cellVolume;
  
  if (cell->IsA("vtkTetra") ){
    cell->GetPoints()->GetPoint(0, r0);
    cell->GetPoints()->GetPoint(1, r1);
    cell->GetPoints()->GetPoint(2, r2);
    cell->GetPoints()->GetPoint(3, r3);
    
    cellVolume = vtkTetra::ComputeVolume(r0, r1, r2, r3);
  }else{
    cout << "ERROR: buggered off to pub before implementing for element type\n" << endl;
    exit(UNKNOWN_CELL_TYPE);
  }
  
  return cellVolume;
}
void CalculateVorticityArray(vtkUnstructuredGrid* grid){
  grid->GetPointData()->SetActiveVectors("Velocity");
  
  vtkCellDerivatives* cellDerivatives = vtkCellDerivatives::New();
  cellDerivatives->SetVectorModeToComputeVorticity();
  cellDerivatives->SetInput(grid);
  cellDerivatives->Update();
  
  grid->GetCellData()->AddArray(cellDerivatives->GetUnstructuredGridOutput()->GetCellData()->GetArray("Vorticity"));
  
  if ( OptionIsSet("debug-vorticity") ){
    grid->GetCellData()->AddArray(cellDerivatives->GetUnstructuredGridOutput()->GetCellData()->GetArray("Tensors"));  
  }
  
  grid->Update();
  cellDerivatives->Delete();
}

void CalculateDebuggingArrays(vtkUnstructuredGrid* grid){
  // Velocity:
  // u = sin(4*pi*y)
  // v = sin(4*pi*x)
  // w = 0
  double xyz[3];
  vtkDataArray* velocity = grid->GetPointData()->GetArray("Velocity");
  double tuple[3];  
  int numberOfCells = grid->GetNumberOfCells();
  vtkDoubleArray* correctVorticity = vtkDoubleArray::New();
  correctVorticity->SetNumberOfComponents(3);
  correctVorticity->SetNumberOfTuples(numberOfCells);  
  
  for (int cellId = 0; cellId < numberOfCells; cellId++){
    // Set up velocity array:
    for (int point = 0; point < 4; point++){
      vtkIdType pointId = grid->GetCell(cellId)->GetPointId(point);      
      grid->GetCell(cellId)->GetPoints()->GetPoint(point, xyz);  
      
      tuple[0] = sin(PIx4 * xyz[1]);// /10.0);
      tuple[1] = sin(PIx4 * xyz[0]);// /10.0);
      tuple[2] = 0;
      velocity->SetTuple(pointId, tuple);    
    }
    
    // Get cell centroid
    vtkCell* cell = grid->GetCell(cellId);
    double parametricCentre[3];
    double cartesianCentre[3];
    int subId = cell->GetParametricCenter(parametricCentre);
    cell->EvaluateLocation(subId, parametricCentre, cartesianCentre, tuple); // Note: vtk reqires 4th parameter, I don't ;)
    
    // Vertical component of vorticity is 4*pi*(cos(4*pi*x) - cos(4*pi*y))
    correctVorticity->SetTuple3( cellId, 0, 0, PIx4 * ( cos(PIx4 * cartesianCentre[0]/10.0) - cos(PIx4 * cartesianCentre[1]/10.0) ) );
  }
  
  correctVorticity->SetName("Known good vorticity");
  grid->GetCellData()->AddArray(correctVorticity);
}


double CalculateIntegral(vtkUnstructuredGrid* grid, Vector3d& result){
  DEBUG("Computing integral...")
    
    // For sanity :)
    result.Zero();
  
  // Each cell should at this point have a vorticity vector.
  vtkDataArray* vorticity = grid->GetCellData()->GetArray("Vorticity");
  
  // Iterate through cells
  double integral=0.0;
  for (vtkIdType i = 0; i < grid->GetNumberOfCells(); i++){
    double cellVolume = fabs(CalculateCellVolume(grid->GetCell(i)));    
    Vector3d cellVorticity( vorticity->GetTuple(i) );
    
    if ( OptionIsSet("2d") ){
      cellVorticity.x = 0;
      cellVorticity.y = 0;
    }
    
    DEBUG("Cell " << i << ": vorticity " << cellVorticity << "  volume: " << cellVolume)
      
      double magnitude = cellVorticity.Magnitude(); 
    integral += magnitude*magnitude*cellVolume;
  }
  return integral;
}

vtkUnstructuredGrid* clip(vtkUnstructuredGrid* grid, string scalar, double value, int clip_above){
  grid->GetPointData()->SetActiveAttribute(scalar.c_str(), vtkDataSetAttributes::SCALARS);
         
  vtkClipDataSet *clipfilter = vtkClipDataSet::New();
  clipfilter->SetInput(grid);
  clipfilter->SetValue(value);
  if(clip_above)
    clipfilter->InsideOutOn();
  else
    clipfilter->InsideOutOff();
  clipfilter->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, scalar.c_str());
  clipfilter->GetOutput()->Update();
  
  return clipfilter->GetOutput();
}

void PrintUsageInformation(const char* executableName){
  cout << "Usage: " << executableName << " -i YourGrid.vtu [OPTIONS]" << endl << endl
    
       << "General options:" << endl
       << "-i, --input=FILENAME" << endl
       << "\tVTU file to use as input" << endl
    
       << "-p, --node-positions" << endl
       << "\tPrint out node XYZ positions" << endl

       << "-b, --buff-body-volume" << endl
       << "\tVolume of buff body" << endl

       << "-c, --clip <scalar name>/<scalar value>/<orientation>" << endl
       << "" << endl

       << "-e, --element-volumes" << endl
       << "\tPrint out element volumes" << endl
    
       << "--vtk-arrays=array1[,array2,...,arrayN]" << endl
       << "\tPrint out contents of specified VTK arrays" << endl
       << "\te.g. --vtk-arrays=Velocity,Temperature" << endl
            
       << "-g, --debug" << endl
       << "\tPrint debugging information" << endl            
    
       << "-o, --offset <scalar name>/<offset value>" << endl
    
       << "-r, --radial-scaling=scale_factor"<<endl
       << "\t Assume shperical earth geometry and scale the radius" << endl
    
       << "-h, --help" << endl
       << "\tPrint this usage information" << endl
            
       << endl
       << "Vorticity integral diagnostic options:" << endl
       
       << "-v, --vorticity-integral" << endl
       << "\tPerform vorticity integral diagnostic" << endl
       
       << "-2, --2d" << endl
       << "\tForce treatment as a 2d problem" << endl       
       
       << "-d, --dump-vtu=FILENAME" << endl
       << "\tDumps a VTU file containing velocity and vorticity to FILENAME" << endl         
       
       << "-w, --debug-vorticity" << endl
       << "\tImposes artificial (sinusoidal) velocity field" << endl
       << "\tand dumps debugging mesh to vorticityDebugMesh.vtu" << endl;       
}

void PrintNodePositions(vtkUnstructuredGrid* grid){
  // Print the positions of all nodes
  cout << "Node positions:" << endl;
  
  double xyz[3];
  for (vtkIdType i = 0; i < grid->GetNumberOfPoints(); i++){
    grid->GetPoint(i, xyz);
    cout << xyz[0] << "\t" << xyz[1] << "\t" << xyz[2] << endl;
  }
  
  cout << "End of node positions" << endl << endl;
}

void TokenizeDelimitedString(const string& stringToTokenize,
           vector<string>& tokenVector,
           const string& delimiters){
  // Sanity check
  if (stringToTokenize.size() == 0) return;
  
  string::size_type leftIndex = 0;
  string::size_type rightIndex = 0;
  
  while (true){
    rightIndex = stringToTokenize.find_first_of(delimiters, leftIndex);
    
    if (rightIndex != string::npos){
      tokenVector.push_back( stringToTokenize.substr(leftIndex, rightIndex - leftIndex) );
      leftIndex = rightIndex + 1;
    }else{
      if (leftIndex != stringToTokenize.size()){
  tokenVector.push_back( stringToTokenize.substr(leftIndex) );
      }
      break;
    }
  }
}

void PrintElementVolumes(vtkUnstructuredGrid* grid){
  for (vtkIdType i=0; i<grid->GetNumberOfCells(); i++)
    cout<<CalculateCellVolume(grid->GetCell(i))<<endl;
}

void PrintVtkArrays(vtkUnstructuredGrid* grid, const string& arrayNameString){
  vector<string> arrayNames;
  TokenizeDelimitedString(arrayNameString, arrayNames, ",");
  
  for (vector<string>::iterator it = arrayNames.begin(); it < arrayNames.end(); it++){
    // Look for point data arrays
    vtkDataArray* array = grid->GetPointData()->GetArray( it->c_str() );
    
    if (array == NULL){
      // Look for cell data arrays
      array = grid->GetCellData()->GetArray( it->c_str() );
    }
    
    if (array != NULL){
      // Print out the array's contents
      cout << "Data from the " << *it << " array:" << endl;
      int numberOfComponents = array->GetNumberOfComponents();
      
      for (int tupleIndex = 0; tupleIndex < array->GetNumberOfTuples(); tupleIndex++){
        for (int componentIndex = 0; componentIndex < numberOfComponents; componentIndex++){
          cout << array->GetComponent(tupleIndex, componentIndex) << "\t";
        }
        
        cout << endl;
      }
      
      cout << "End of data from the " << *it << " array" << endl;
    }else{
      cout << "Unable to find array \"" << *it << "\" in this VTU file" << endl;
    }
  }
}

void ParseArguments(int argc, char** argv){
  struct option longOptions[] = {
    {"2d", 0, 0, '2'},
    {"buff-body-volume", 0, 0, 'b'},
    {"clip", 0, 0, 'c'},
    {"debug", 0, 0, 'g'},
    {"debug-vorticity", 0, 0, 'w'},
    {"dump-vtu", 1, 0, 'd'},
    {"element-volumes", 0, 0, 'e'},
    {"help", 0, 0, 'h'},
    {"input", 1, 0, 'i'},
    {"node-positions", 0, 0, 'p'},
    {"offset", 0, 0, 'o'},
    {"radial-scaling", 1, 0, 'r'},
    {"vorticity-integral", 0, 0, 'v'},
    {"vtk-arrays", 1, 0, 'a'},
    {0, 0, 0, 0}
  };
  
  int optionIndex = 0;
  int c;
  
  while (true){
    c = getopt_long(argc, argv, "bc:i:2vd:r:hgwpeao:", longOptions, &optionIndex);
    if (c == -1) break;
    
    switch (c){  
    case '2':
      optionMap["2d"] = "true";
      break;        

    case 'a':
      optionMap["vtk-arrays"] = optarg;
      break;

    case 'b':
      optionMap["buff-body-volume"] = "true";
      break;

    case 'c':
      optionMap["clip"] = optarg;
      break;

    case 'd':
      optionMap["dump-vtu"] = optarg;
      break;

    case 'e':
      optionMap["element-volumes"] = "true";
      break;  

    case 'g':
      optionMap["debug"] = "true";
      break;

    case 'h':
      PrintUsageInformation(argv[0]);
      exit(SUCCESS);

    case 'i':
      optionMap["input"] = optarg;
      break;
            
    case 'o':
      optionMap["offset"] = optarg;
      break;

    case 'p':
      optionMap["node-positions"] = "true";
      break;

    case 'r':
      optionMap["radial-scaling"] = optarg;
      break;

    case 'v':
      optionMap["vorticity-integral"] = "true";
      break;  

    case 'w':
      optionMap["debug-vorticity"] = "true";
      break;
      
    case '?':
      PrintUsageInformation(argv[0]);
      exit(BAD_ARGUMENTS);
    }
  }
  
  if (optionMap.find("input") == optionMap.end()){
    cout << "No input file specified!" << endl;
    exit(BAD_ARGUMENTS);
  }
}

double buff_body_volume(vtkUnstructuredGrid* grid){
  double volume = 0.0;
  for (vtkIdType i=0; i<grid->GetNumberOfCells(); i++)
    volume += fabs(CalculateCellVolume(grid->GetCell(i)));

  double bbox[6];
  grid->GetBounds(bbox);
  double bbox_volume = (bbox[1]-bbox[0])*(bbox[3]-bbox[2])*(bbox[5]-bbox[4]);
  
  return bbox_volume - volume;
}

int main(int argc, char** argv){
  ParseArguments(argc, argv);
  
  // Load
  vtkUnstructuredGrid* grid = LoadGrid( optionMap["input"] );

  bool wrote_vtu=false;

  // Buff body volume?
  if (OptionIsSet("buff-body-volume")){
    cout<<"Volume of buff body (i.e. void region): "<<buff_body_volume(grid)<<endl;
  }

  // Clip mesh
  if (OptionIsSet("clip")){
    vtkUnstructuredGrid* clipped_data;
    vector<string> tokens;
    TokenizeDelimitedString(optionMap["clip"], tokens, "/");
    for(size_t i=0;i<tokens.size();i+=3){
      if((i+2)==tokens.size()){
        cerr<<"No upper limit specified for clip.\n";
        PrintUsageInformation(argv[0]);
        exit(BAD_ARGUMENTS);
      }
      double svalue = strtod(tokens[i+1].c_str(), NULL);
      int clip_above = strtol(tokens[i+2].c_str(), NULL, 0);
      clipped_data = clip(grid, tokens[i], svalue, clip_above);
    }
    
    WriteGrid(clipped_data, "clipped.vtu");
  }
  
  // Offset scalar value
  if (OptionIsSet("offset")){
    vtkUnstructuredGrid* offset_data = vtkUnstructuredGrid::New();
    vector<string> tokens;
    TokenizeDelimitedString(optionMap["offset"], tokens, "/");
    if(tokens.size()!=2){
      cerr<<"Bad arguments for offset option.\n";
      PrintUsageInformation(argv[0]);
      exit(BAD_ARGUMENTS);
    }
    double offset = strtod(tokens[1].c_str(), NULL);

    offset_data->SetPoints(grid->GetPoints());
    offset_data->SetCells(grid->GetCellType(0), grid->GetCells());
    offset_data->GetPointData()->SetScalars(grid->GetPointData()->GetScalars(tokens[0].c_str()));
    
    for(int i=0;i<offset_data->GetNumberOfPoints(); i++){
      double nscalar = offset_data->GetPointData()->GetScalars(tokens[0].c_str())->GetTuple1(i)+offset;
      offset_data->GetPointData()->GetScalars(tokens[0].c_str())->SetTuple1(i, nscalar);
    }

    WriteGrid(offset_data, "offset.vtu");
  }
  
  // Process
  if ( OptionIsSet("vorticity-integral") ){
    if ( OptionIsSet("debug-vorticity") ){
      CalculateDebuggingArrays(grid);
      CalculateVorticityArray(grid);
      WriteGrid(grid, "vorticityDebugMesh.vtu");
      wrote_vtu=true;
    }
    
    CalculateVorticityArray(grid);
    
    if ( OptionIsSet("dump-vtu") ){
      WriteGrid( grid, optionMap["dump-vtu"] );
      wrote_vtu=true;
    }
    
    // Find integral
    Vector3d result;
    double integral = CalculateIntegral(grid, result);
    
    cout << integral << endl;    
  }else{
    if ( OptionIsSet("debug-vorticity") || OptionIsSet("dump-vtu") ){
      cout << "Please use the --vorticity-integral (-v) option to "
     << "enable the vorticity diagnostic" << endl;
      return BAD_ARGUMENTS;
    }
  }
  
  if ( OptionIsSet("node-positions") ){
    PrintNodePositions(grid);
  }
  
  if ( OptionIsSet("vtk-arrays") ){
    PrintVtkArrays(grid, optionMap["vtk-arrays"]);
  }

  if( OptionIsSet("element-volumes") ){
    PrintElementVolumes(grid);
  }
  
  if((!wrote_vtu)&&OptionIsSet("radial-scaling")){
    WriteGrid(grid, "radial.vtu");
    wrote_vtu=true;
  }

  if(!wrote_vtu){
    WriteGrid(grid, "dump.vtu");
    wrote_vtu=true;
  }

  return SUCCESS;
}
