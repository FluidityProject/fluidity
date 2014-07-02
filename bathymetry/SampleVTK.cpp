#include <cassert>

#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <algorithm>

#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkCell.h>
#include <vtkPointLocator.h>
#include <vtkPolyData.h>
#include <vtkProbeFilter.h>
#include <vtkIdTypeArray.h>
#include <vtkPointData.h>

extern "C" {
int sample_vtu(const char *filename, const char *fieldname, const double *x, const double *y, int cnt,
	       double *values){
  
  // Read in file.
  vtkXMLUnstructuredGridReader *reader = vtkXMLUnstructuredGridReader::New();
  reader->SetFileName(filename);
  reader->Update();
  
  vtkUnstructuredGrid *ug = reader->GetOutput();
  
  // Initialise point locator
  vtkPointLocator *locator = vtkPointLocator::New();
  locator->SetDataSet(ug);
  locator->SetTolerance(10.0);
  locator->Update();
  
  // Initialise probe
  vtkPoints *points = vtkPoints::New();
  for(int i=0;i<cnt;i++)
    points->InsertNextPoint(x[i], y[i], 0);
  
  vtkPolyData *polydata = vtkPolyData::New();
  polydata->SetPoints(points);
  
  vtkProbeFilter *probe = vtkProbeFilter::New();
  probe->SetInput(polydata);
  probe->SetSource(ug);
  probe->Update();
  
  // Get the values
  vtkIdTypeArray *valid_ids = probe->GetValidPoints();
  vtkPointData *pd = probe->GetOutput()->GetPointData();
  int valid_loc = 0;
  int invalid_cnt = 0;
  for(int i=0;i<cnt;i++){
    if(valid_ids->GetTuple1(valid_loc) == i){
      values[i] = pd->GetArray(fieldname)->GetTuple1(i);
      valid_loc++;
    }else{
      invalid_cnt++;
      double coordinate[] = {x[i], y[i], 0.0};
      int nearest = locator->FindClosestPoint(coordinate);
      values[i] = ug->GetPointData()->GetArray(fieldname)->GetTuple1(nearest);
    }
  }

  return invalid_cnt;
}
}
// Compile using  g++ -I/usr/include/vtk-5.8 -Wno-deprecated SampleVTK.cpp -lvtkCommon  -lvtkIO -lvtkFiltering -lvtkGenericFiltering -lvtkGraphics

// int main(){
//   std::string filename("bathymetry.vtu");
  
//   // We are going to take the X,Y values for the probe from the input
//   // file to show that we can find the points exactly.
//   vtkXMLUnstructuredGridReader *reader = vtkXMLUnstructuredGridReader::New();
//   reader->SetFileName(filename.c_str());
//   reader->Update();
  
//   vtkUnstructuredGrid *ug = reader->GetOutput();
//   size_t NPoints = ug->GetNumberOfPoints();
//   std::vector<double> x(NPoints), y(NPoints);
  
//   for(size_t i=0;i<NPoints;i++){
//     const double *point = ug->GetPoint(i);
//     x[i] = point[0];
//     y[i] = point[1];
//   }

//   // Test
//   std::vector<double> values(NPoints); 
//   int invalid_cnt = sample_vtu(filename.c_str(), "BOTTOM", &(x[0]), &(y[0]), NPoints, &(values[0]));
  
//   std::cout<<"Interpolated "<<NPoints<<" values. Of these "
// 	   <<invalid_cnt<<" are nearest-neighbour points."<<std::endl;

//   std::cout<<"First 10 values"<<std::endl;
//   for(size_t i=0;i<std::min((size_t)10, NPoints);i++){
//     double s = ug->GetPointData()->GetArray("BOTTOM")->GetTuple1(i);
//     std::cout<<"source="<<s<<", interpolated="<<values[i]<<", diff="<<s-values[i]<<std::endl;
//   }

//   reader->Delete();
//   return 0;
// }
