/*  Copyright (C) 2009 Imperial College London and others.
    
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
#include "confdefs.h"
#include "NEMOReader.h"

using namespace std;

NEMOReader NEMOReader_v2_global;

extern int projections(int nPoints, double *x, double *y, double *z, string current_coord, string output_coord);


extern "C" {
#define rotate_ll2cart_fc F77_FUNC(rotate_ll2cart, ROTATE_LL2CART)
  extern void rotate_ll2cart_fc(double *longitude, double *latitude, double *u, double *v,
                     const double *r3u, const double *r3v, const double *r3w);

#define get_nemo_variables_fc F77_FUNC(get_nemo_variables, GET_NEMO_VARIABLES)
    void get_nemo_variables_fc(double *time, const double *X, const double *Y, const double *Z, 
                     double *Te, double *Sa, double *U, double *V, double *W, double *SSH, const int *NNodes);
}

void get_nemo_variables_fc(double *time, const double *X, const double *Y, const double *Z, 
                     double *Te, double *Sa, double *U, double *V, double *W, double *SSH, const int *n){

    const int nFields = 5; // number of fields being read in, currently temperature, salinity, U, V and sea surface height
    const int NNodes = *n;
        
    double *salinity = new double[NNodes];
    double *seatemp = new double[NNodes];
    double *uvel = new double[NNodes];
    double *vvel = new double[NNodes];
    double *wvel = new double[NNodes];
    double *seash = new double[NNodes];

    NEMOReader_v2_global.SetTimeSeconds(*time);
    
    // convert from Cart to long-lat
    double *x = new double[NNodes];
    double *y = new double[NNodes];
    double *z = new double[NNodes];
    double *depth = new double[NNodes];
    for (int i = 0; i < NNodes; i++) {
        x[i] = X[i];
        y[i] = Y[i];
        z[i] = Z[i];
        depth[i]=6378000.0-sqrt(X[i]*X[i]+Y[i]*Y[i]+Z[i]*Z[i]);
    }    

    int ret = projections(NNodes, x, y, z, "cart", "spherical");
    if (ret != 0) {
        cerr<<"Error converting coord system"<<endl;
    }
    
    /*
     * The "values" array below contains the scalar
     * values from the NEMOdata netCDF file. The table
     * below shows which NEMO parameter is at which
     * array index, along with a physical meaning
     *
    NEMO field  | Index  | Physical meaning   |     units       |
    ------------+-------+---------------------+-----------------+
    temperature |   0   | Sea Temperature     | degrees celcius |
    salinity    |   1   | Salinity            | PPT             |
    u           |   2   | azimuthal velocity  | ms^-1           |
    v           |   3   | meridional velocity | ms^-1           |
    ssh         |   4   | Sea surface height  | m               |
    */

    // loop over nodes
    for (int i=0; i<NNodes; i++) {
        
        double latitude = y[i]; 
        double longitude = x[i];
        double depth_p = depth[i];

        // get values from the netcdf file
        double values[nFields];

        NEMOReader_v2_global.GetScalars(longitude, latitude, depth_p, values);

        // values contains the values above in the order we registered them
        // See above
        seatemp[i] = values[0];
        salinity[i] = values[1];
        seash[i] = values[4];

        // rotate wind to cartesian grid
        double u_rot = values[2];
        double v_rot = values[3];

        // This currently assumes that you are running on the sphere. If NEMO data is going to be used for
        // planar cases an if if(have_option("/geometry/spherical_earth")) will need to be passed or
        // added.
        rotate_ll2cart_fc(&longitude, &latitude, &u_rot, &v_rot, &uvel[i], &vvel[i], &wvel[i]);
        
    }

    delete [] x;
    delete [] y;
    delete [] z;
    delete [] depth;

    for (int i=0; i<NNodes; i++) {

      Te[i]=seatemp[i];
      Sa[i]=salinity[i];
      U[i]=uvel[i];
      V[i]=vvel[i];
      W[i]=wvel[i];
      SSH[i]=seash[i];

    }

  // clean up
  delete [] salinity;
  delete [] seatemp;
  delete [] uvel;
  delete [] vvel;
  delete [] wvel;
  delete [] seash;

}

#ifdef NEMODATALOAD_UNITTEST
#include <vector>
#include <vtk.h>

void vtk_add_scalar(vector<double> &scalar, const char *scalar_name, vtkUnstructuredGrid *ug){
  int NNodes = ug->GetNumberOfPoints();
  
  vtkDoubleArray *newScalars = vtkDoubleArray::New();
  newScalars->SetName(scalar_name);
  newScalars->SetNumberOfComponents(1);
  newScalars->SetNumberOfTuples(NNodes);
  
  for(int i=0; i<NNodes; i++)
    newScalars->InsertValue(i, scalar[i]);
  
  ug->GetPointData()->AddArray(newScalars);
  ug->GetPointData()->SetActiveAttribute(scalar_name, vtkDataSetAttributes::SCALARS);
  
  ug->Update();
  newScalars->Delete();
}

void vtk_add_vector(vector<double> &vx, vector<double> &vy, vector<double> &vz, const char *vector_name, vtkUnstructuredGrid *ug){
  int NNodes = ug->GetNumberOfPoints();
  
  vtkDoubleArray *newVectors = vtkDoubleArray::New();
  newVectors->SetName(vector_name);
  newVectors->SetNumberOfComponents(3);
  newVectors->SetNumberOfTuples(NNodes);
  
  for(int i=0; i<NNodes; i++)
    newVectors->SetTuple3(i, vx[i], vy[i], vz[i]);;
  
  ug->GetPointData()->AddArray(newVectors);
  ug->GetPointData()->SetActiveAttribute(vector_name, vtkDataSetAttributes::VECTORS);
  
  ug->Update();
  newVectors->Delete();
}

// This is a simple test program to view the interpolation of the NEMO data onto an
// unstructured grid supplied from a vtu.
int main(int argc, char **argv){

    char *filename=argv[1];
    vtkXMLUnstructuredGridReader *reader = vtkXMLUnstructuredGridReader::New();
    reader->SetFileName(filename);
    reader->Update();
  
    vtkUnstructuredGrid *ug = reader->GetOutput();

    int NNodes = ug->GetNumberOfPoints();
    double time = 0.0;  
    vector<double> Te(NNodes, 0.0), Sa(NNodes, 0.0), U(NNodes, 0.0), V(NNodes, 0.0), 
    W(NNodes, 0.0), SSH(NNodes, 0.0), X(NNodes, 0.0), Y(NNodes, 0.0), Z(NNodes, 0.0), 
    x(NNodes,0.0), y(NNodes,0.0), z(NNodes,0.0);
 
    for (int i=0; i<NNodes; i++) {
        double r[3];
        ug->GetPoints()->GetPoint(i, r);
        X[i] = r[0];
        Y[i] = r[1];
        Z[i] = r[2];
        // these will be turned into lat/long
        x[i] = X[i];
        y[i] = Y[i];
        z[i] = Z[i];
    }

    NEMOReader NEMOReader_NEMOdata;

    // this data set is the input
    NEMOReader_v2_global.RegisterDataFile("/data/ORCA_grid/gmt_gridding/NEMOdata.nc");

    NEMOReader_v2_global.SetSimulationTimeUnits("seconds since 1987-01-05 00:00:0.0");
    NEMOReader_v2_global.AddFieldOfInterest("temperature");  //  0   | Sea temperature
    NEMOReader_v2_global.AddFieldOfInterest("salinity");     //  1   | Salinity
    NEMOReader_v2_global.AddFieldOfInterest("u");            //  2   | Azimuthal velocity
    NEMOReader_v2_global.AddFieldOfInterest("v");            //  3   | Meridional velocity
    NEMOReader_v2_global.AddFieldOfInterest("ssh");          //  4   | Sea surface height
    
  
    int ret = projections(NNodes, &x[0], &y[0], &z[0], "cart", "spherical");
    if (ret != 0) {
        cerr<<"Error converting coord system"<<endl;
    }
    int const nFields_in=4;
    double values_in[nFields_in];

    // now calucalte fluxes for the whole globe and put on the output VTU file for
    // visualisation

    // Note that latitudinal velocities are currently missing from the NEMO data
    // and is thus set to zero.

    get_nemo_variables_fc(&time, &X[0], &Y[0], &Z[0], &Te[0], &Sa[0], &U[0], &V[0], &W[0],
                        &SSH[0], &NNodes);

    // add to VTU file
    vtk_add_scalar(Te, "Temperature", ug);
    vtk_add_scalar(Sa, "Salinity", ug);
    vtk_add_scalar(U, "U", ug);
    vtk_add_scalar(V, "V", ug);
    vtk_add_scalar(W, "W", ug);
    vtk_add_vector(U, V, W, "Velocity", ug);
    vtk_add_scalar(SSH, "Sea surface height", ug);
    
    // now print values at a single node
    int n = 1;
    // A point in the North West Approaches
    double londum,latdum,REdum,pidum;
    pidum=acos(-1.0);
    REdum=6378100.0;
    londum=(-8.0)*pidum/180.0;
    latdum=(60.0)*pidum/180.0;
    X[0] = REdum*cos(latdum)*cos(londum);
    Y[0] = REdum*cos(latdum)*sin(londum);
    Z[0] = REdum*sin(latdum);
    double location_depth=6378100.0-sqrt(X[0]*X[0]+Y[0]*Y[0]+Z[0]*Z[0]);
    y[0] = Y[0] , x[0] = X[0] , z[0] = Z[0];
    ret = projections(n, &x[0], &y[0], &z[0], "cart", "spherical");
    if (ret != 0) {
        cerr<<"Error converting coord system"<<endl;
    }
    
    NEMOReader_v2_global.GetScalars(x[0], y[0], location_depth, values_in);

    Te[0] = values_in[0];
    Sa[0] = values_in[1];
    SSH[0] = values_in[4];

    // rotate wind to cartesian grid
    double u_sphere = values_in[2];
    double v_sphere = values_in[3];

    // This currently assumes that you are running on the sphere. If NEMO data is going to be used for
    // planar cases an if if(have_option("/geometry/spherical_earth")) will need to be passed or
    // added.
    rotate_ll2cart_fc(&londum, &latdum, &u_sphere, &v_sphere, &U[0], &V[0], &W[0]);

    cout <<endl<<"Values for "<<x[0]<<","<<y[0]<<endl<<endl;

    cout <<"    Variable         |     NEMO value     |"<<endl;
    cout <<"---------------------+--------------------|"<<endl;
    cout <<"    Temperature      |    "<<Te[0] << endl;
    cout <<"     Salinity        |    "<<Sa[0] << endl;
    cout <<"        U            |    "<<U[0] << endl;
    cout <<"        V            |    "<<V[0] << endl;
    cout <<"        W            |    "<<W[0] << endl;
    cout <<"       SSH           |    "<<SSH[0] << endl;

  vtkXMLUnstructuredGridWriter *writer = vtkXMLUnstructuredGridWriter::New();
  writer->SetFileName("NEMOvalues.vtu");
  writer->SetInput(ug);
  writer->Write();
}
#endif