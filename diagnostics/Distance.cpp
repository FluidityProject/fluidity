#include "Distance.h"

#include "vtkCell.h"
#include "vtkLine.h"
#include "vtkPlane.h"
#include "vtkGenericCell.h"
#include "vtkCellLocator.h"
#include "vtkSmartPointer.h"

#include <math.h>
#include <limits>
#include <algorithm>

Distance* Distance::New(){ return new Distance;}
void Distance::Delete(){ delete this;}

Distance::Distance(){this->boundary=NULL;};
Distance::~Distance(){}

void Distance::SetBoundary(vtkUnstructuredGrid* b) {
  if(this->boundary) this->boundary->Delete();
  this->boundary = b;
}

vtkUnstructuredGrid* Distance::GetBoundary(){
  return this->boundary;
}

int Distance::CalculateDistance(vtkUnstructuredGrid* mesh,
				vtkUnstructuredGrid* boundary,
				double* distance) {

  vtkCellLocator *locator= vtkCellLocator::New();
  locator->SetDataSet(boundary);
  locator->BuildLocator();
  double x[3], cp[3];
  vtkIdType cellId;
  int subId;

  for (vtkIdType i=0; i<mesh->GetNumberOfPoints(); ++i) {

    mesh->GetPoint(i,x);
    vtkSmartPointer<vtkGenericCell> c= vtkGenericCell::New();
    locator->FindClosestPoint(x, cp, c, cellId, subId, distance[i]);

    distance[i] = sqrt(distance[i]);
    c->Delete();
  }

  locator->Delete();

  return 1;
}

int Distance::CalculateDirection(vtkUnstructuredGrid* mesh,
				vtkUnstructuredGrid* boundary,
				int dim, double* direction) {

  vtkCellLocator *locator= vtkCellLocator::New();
  locator->SetDataSet(boundary);
  locator->BuildLocator();
  double x[3], cp[3];
  vtkIdType cellId;
  int subId;

  for (vtkIdType i=0; i<mesh->GetNumberOfPoints(); ++i) {

    double dist = std::numeric_limits<double>::infinity();
    mesh->GetPoint(i,x);

    vtkSmartPointer<vtkGenericCell> c= vtkGenericCell::New();
    locator->FindClosestPoint(x, cp, c, cellId, subId, dist);

    for (int j=0;j<dim;++j) {
      direction[dim*i+j] = x[j]-cp[j];
    }
    c->Delete();
    
  }

  locator->Delete();

  return 1;
}

extern "C" {

  void vtk_distance(int dim, int nodes, int snodes, int eles, int seles,
		    double* pts, int* ndglno, double* spts, int* sndglno,
		    double* distance) {
    vtkPoints* vtkpts = vtkPoints::New();
    vtkpts->SetNumberOfPoints(nodes);
    double x[3] ={0.,0.,0.};
    for (vtkIdType i=0;i<nodes;++i){
      for(int j=0;j<dim;++j){
	x[j] = pts[dim*i+j];
      }
      vtkpts->SetPoint(i,x);
    }
    vtkPoints* svtkpts = vtkPoints::New();
    svtkpts->SetNumberOfPoints(snodes);
    for (vtkIdType i=0;i<snodes;++i){
      for(int j=0;j<dim;++j){
	x[j] = spts[dim*i+j];
      }
      svtkpts->SetPoint(i,x);
    }
    vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::New();
    vtkUnstructuredGrid* boundary = vtkUnstructuredGrid::New();
    ugrid->SetPoints(vtkpts);
    vtkIdList* list = vtkIdList::New();
    list->SetNumberOfIds(dim+1);
    int type[4] = {VTK_PIXEL, VTK_LINE, VTK_TRIANGLE, VTK_TETRA};
    for (vtkIdType i=0; i<eles; ++i) {
      for(int j=0; j<dim+1;++j) {
	list->SetId(j, ndglno[(dim+1)*i+j]-1);
      }
      ugrid->InsertNextCell(type[dim], list);
    }
    boundary->SetPoints(svtkpts);
    vtkIdList* list2 = vtkIdList::New();
    list2->SetNumberOfIds(dim);
    for (vtkIdType i=0; i<seles; ++i) {
      for(int j=0; j<dim;++j) {
	list2->SetId(j, sndglno[(dim)*i+j]-1);
      }
      boundary->InsertNextCell(type[dim-1], list2);
    }
    Distance* df = Distance::New();
    df->CalculateDistance(ugrid, boundary, distance);

    df->Delete();
    boundary->Delete();
    ugrid->Delete();
    vtkpts->Delete();

  }


  void vtk_direction(int dim, int nodes, int snodes, int eles, int seles,
		    double* pts, int* ndglno, double* spts, int* sndglno,
		    double* direction) {
    vtkPoints* vtkpts = vtkPoints::New();
    vtkpts->SetNumberOfPoints(nodes);
    double x[3] ={0.,0.,0.};
    for (vtkIdType i=0;i<nodes;++i){
      for(int j=0;j<dim;++j){
	x[j] = pts[dim*i+j];
      }
      vtkpts->SetPoint(i,x);
    }
    vtkPoints* svtkpts = vtkPoints::New();
    svtkpts->SetNumberOfPoints(snodes);
    for (vtkIdType i=0;i<snodes;++i){
      for(int j=0;j<dim;++j){
	x[j] = spts[dim*i+j];
      }
      svtkpts->SetPoint(i,x);
    }
    vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::New();
    vtkUnstructuredGrid* boundary = vtkUnstructuredGrid::New();
    ugrid->SetPoints(vtkpts);
    vtkIdList* list = vtkIdList::New();
    list->SetNumberOfIds(dim+1);
    int type[4] = {VTK_PIXEL, VTK_LINE, VTK_TRIANGLE, VTK_TETRA};
    for (vtkIdType i=0; i<eles; ++i) {
      for(int j=0; j<dim+1;++j) {
	list->SetId(j, ndglno[(dim+1)*i+j]-1);
      }
      ugrid->InsertNextCell(type[dim], list);
    }
    boundary->SetPoints(svtkpts);
    vtkIdList* list2 = vtkIdList::New();
    list2->SetNumberOfIds(dim);
    for (vtkIdType i=0; i<seles; ++i) {
      for(int j=0; j<dim;++j) {
	list2->SetId(j, sndglno[(dim)*i+j]-1);
      }
      boundary->InsertNextCell(type[dim-1], list2);
    }
    Distance* df = Distance::New();
    df->CalculateDirection(ugrid, boundary, dim, direction);

    df->Delete();
    boundary->Delete();
    ugrid->Delete();
    vtkpts->Delete();
    svtkpts->Delete();

  }

}
