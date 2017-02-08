//    Copyright (C) 2006 Imperial College London and others.
//   
//    Please see the AUTHORS file in the main source directory for a full list
//    of copyright holders.
//
//    Prof. C Pain
//    Applied Modelling and Computation Group
//    Department of Earth Science and Engineering
//    Imperial College London
//
//    amcgsoftware@imperial.ac.uk
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation,
//    version 2.1 of the License.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
//    USA

#include "vtkCell.h"
#include "vtkDataArray.h"
#include "vtkPoints.h"
#include "vtkCellData.h"
#include "vtkDataSet.h"
#include "vtkUnstructuredGrid.h"
#include "vtkMeshQuality.h"

#include <stdio.h>


extern "C" {
void mesh_quality_c(int* dim, int* n_nodes, int* n_elements, int* connectivity_len,
		    int* measure, double* points, int* connectivity, double* quality) {
  vtkPoints* pts = vtkPoints::New();
  vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::New();
  pts->SetDataTypeToDouble();
  pts->Allocate(*n_nodes);
  double x[3]={0.0,0.0,0.0};
  for (vtkIdType i=0; i<*n_nodes; ++i) {
    for (int j=0;j<*dim; j++) {
      x[j] = points[(*dim)*i+j];
    }
    pts->SetPoint(i,x);
  }
  ugrid->SetPoints(pts);

  int cell_type;
  vtkIdType cell[20];
  if (*dim==2) {
    cell_type = VTK_TRIANGLE;
  } else if (*dim==3) {
    cell_type = VTK_TETRA;
  }

  ugrid->Allocate(0);

  for(vtkIdType i=0;i<*n_elements;++i) {
    for (int j=0;j<*dim+1;++j) {
      // off by one due to fortran indexing
      cell[j] = connectivity[(*dim+1)*i+j]-1;
    }
    ugrid->InsertNextCell(cell_type,*dim+1,cell);
  }

  vtkMeshQuality* filter = vtkMeshQuality::New();

#if VTK_MAJOR_VERSION <= 5
  filter->SetInput(ugrid);
#else
  filter->SetInputData(ugrid);
#endif
  
  filter->SetTriangleQualityMeasure(*measure);
  if (*measure == VTK_QUALITY_AREA) {
    filter->SetTetQualityMeasure(VTK_QUALITY_VOLUME);
  } else {
    filter->SetTetQualityMeasure(*measure);
  }
  filter->Update();
  vtkDataSet* vgrid=filter->GetOutput();
  vtkDataArray* data = vgrid->GetCellData()->GetArray(0);

  for (vtkIdType i=0;i<data->GetNumberOfTuples();++i) {
    quality[i] = data->GetTuple1(i);
  }

  filter->Delete();
  ugrid->Delete();
  pts->Delete();

}

}
