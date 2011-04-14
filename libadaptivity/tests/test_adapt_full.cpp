/* Copyright (C) 2009 Imperial College London.

 Please see the AUTHORS file in the main source directory for a full
 list of copyright holders.

 Dr Gerard J Gorman
 Applied Modelling and Computation Group
 Department of Earth Science and Engineering
 Imperial College London

 g.gorman@imperial.ac.uk

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
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
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "vtk.h"

#include "ErrorMeasure.h"
#include "Adaptivity.h"
#include "DiscreteGeometryConstraints.h"

/** This test seeks to execute a full adaptive mesh example.
 */

using namespace std;

int main(int argc, char **argv){
  vtkXMLUnstructuredGridReader *ug_reader = vtkXMLUnstructuredGridReader::New();
  ug_reader->SetFileName("output.vtu");
  ug_reader->Update();
  
  vtkUnstructuredGrid *ug = ug_reader->GetOutput();
  ug->Update();
  
  // These should be populated only once at the outset of a simulation
  // and be maintained thereafter.
  vector<int> SENList, sids;
  {
    DiscreteGeometryConstraints constraints;
    constraints.verbose_on();
    constraints.set_coplanar_tolerance(0.9999999);
    constraints.set_volume_input(ug);
    
    constraints.get_coplanar_ids(sids);
    constraints.get_surface(SENList);
    assert(sids.size()*3==SENList.size());
    cout<<"Found "<<sids.size()<<" surface elements\n";
  }

  for(size_t i=0;i<1;i++){
    cout<<"Outer iteration: "<<i<<endl;
    
    DiscreteGeometryConstraints constraints;
    constraints.verbose_on();
    constraints.set_surface_input(ug, SENList, sids);
    
    vector<double> max_len;
    constraints.get_constraints(max_len);
    constraints.write_vtk(string("sids.vtu"));
    
    /* Test merging of metrics.
     */
    ErrorMeasure error;
    error.verbose_on();
    error.set_input(ug);
    error.add_field("Vm", 1.0, false, 0.01);
    error.set_max_length(2.0);
    error.set_max_length(&(max_len[0]), ug->GetNumberOfPoints());
    error.set_min_length(0.002);
    error.apply_gradation(1.3);
    error.set_max_nodes(200000);
    
    error.diagnostics();
    
    vtkXMLUnstructuredGridWriter *metric_writer = vtkXMLUnstructuredGridWriter::New();
    metric_writer->SetFileName("metric.vtu");
    metric_writer->SetInput(ug);
    metric_writer->Write();
    metric_writer->Delete();
    
    ug->GetPointData()->RemoveArray("mean_desired_lengths");
    ug->GetPointData()->RemoveArray("desired_lengths");
    
    Adaptivity adapt;
    adapt.verbose_on();
    adapt.set_from_vtk(ug, true);
    adapt.set_adapt_sweeps(5);
    adapt.set_surface_mesh(SENList);
    adapt.set_surface_ids(sids);
    adapt.adapt();
    adapt.get_surface_ids(sids);
    adapt.get_surface_mesh(SENList);
    vtkUnstructuredGrid *tmp_ug = adapt.get_adapted_vtu();
    if(i>0)
      ug->Delete();
    ug = tmp_ug;

    ug->GetPointData()->RemoveArray("metric");

  }
  
  vtkXMLUnstructuredGridWriter *ug_writer = vtkXMLUnstructuredGridWriter::New();
  ug_writer->SetFileName("adapted.vtu");
  ug_writer->SetInput(ug);
  ug_writer->Write();

  ug_writer->Delete();
  ug_reader->Delete();
  return 0;
}
