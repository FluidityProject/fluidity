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
#include "FluxesReader.h"

using namespace std;

//FluxesReader FluxesReader_global;

extern int projections(int nPoints, double *x, double *y, double *z, string current_coord, string output_coord);


extern "C" {
#define ncar_ocean_fluxes_fc F77_FUNC(ncar_ocean_fluxes, NCAR_OCEAN_FLUXES)
  extern void ncar_ocean_fluxes_fc(const int *NNodes, double *speed, double *t, double *SST,
                                   double *q, double *qs, const double *z, bool *avail,
                                   double *cd, double *ch, double *ce, 
                                   double *ustar, double *bstar);

#define rotate_wind_fc F77_FUNC(rotate_wind, ROTATE_WIND)
  extern void rotate_wind_fc(double *longitude, double *latitude, const double *r3u, const double *r3v, const double *r3w,
                            double *u, double *v);

#define get_era40_fluxes_fc F77_FUNC(get_era40_fluxes, GET_ERA40_FLUXES)
    void get_era40_fluxes_fc(double *time, const double *X, const double *Y, const double *Z, 
                     double *T, const double *Vx, const double *Vy, const double *Vz,
                     double *F_as, double *Q_as, double *tau_u, double *tau_v, double *Q_solar,
                     const int *NNodes, bool rotate, double *q_l = NULL, double *q_h = NULL, double *e = NULL,
                     double *q_e = NULL, double *q_p = NULL);
}

void get_era40_fluxes_fc(double *time, const double *X, const double *Y, const double *Z, 
                     double *T, const double *Vx, const double *Vy, const double *Vz,
                     double *F_as, double *Q_as, double *tau_u, double *tau_v, double *Q_solar,
                     const int *n, bool rotate, double *q_l, double *q_h, double *e, double *q_e, double *q_p ){

    // physical constants
    const double air_density = 1.22, 
               vapour_latent = 2.5e6, 
               air_specificHeat = 1000.5, 
               alpha = 0.066, 
               sb = 5.67e-8, 
               fusion_latent = 3.337e5, 
               q1 = 0.98*640380, 
               q2 = -5107.4,
               ocean_density = 1027.0,
               ocean_heat_capacity = 4000.0,
               ref_salinity = 35.0,
               accumulated_correction = 6.0*60.0*60.0;

    // problem constants
    const int nFields = 9;
    const int NNodes = *n;
        
    double *delU_u = new double[NNodes];
    double *delU_v = new double[NNodes];
    double *cd = new double[NNodes];
    double *ch = new double[NNodes];
    double *ce = new double[NNodes];
    double *t_2m = new double[NNodes];
    double *q = new double[NNodes];
    double *z_data = new double[NNodes];
    double *solar = new double[NNodes]; 
    double *thermal = new double[NNodes];
    double *cloud = new double[NNodes];
    double *runoff = new double[NNodes];
    double *ppt = new double[NNodes];
    double *speed = new double[NNodes];
    double *qs = new double[NNodes];
    double *bstar = new double[NNodes];
    double *ustar = new double[NNodes];
    double *SST = new double[NNodes];
    bool *avail = new bool[NNodes];

    FluxesReader_global.SetTimeSeconds(*time);
    
    // convert from Cart to long-lat
    double *x = new double[NNodes];
    double *y = new double[NNodes];
    double *z = new double[NNodes];
    for (int i = 0; i < NNodes; i++) {
        x[i] = X[i];
        y[i] = Y[i];
        z[i] = Z[i];
    }

    int ret = projections(NNodes, x, y, z, "cart", "spherical");
    if (ret != 0) {
        cerr<<"Error converting coord system"<<endl;
    }
    
    /*
     * The "values" array below contains the scalar
     * values from the ERA40 netCDF file. The table
     * below shows which ERA40 parameter is at which
     * array index, along with a physical meaning
     *
    ERA40 field | Index | Physical meaning
    ------------+-------+-----------------
    10u         |   0   | 10 metre U wind component
    10v         |   1   | 10 metre V wind component
    ssrd        |   2   | Surface solar radiation
    strd        |   3   | Surface thermal radiation 
    ro          |   4   | Runoff
    tp          |   5   | Total precipitation
    2d          |   6   | Dewpoint temp at 2m
    2t          |   7   | Air temp at 2m
    msl         |   8   | Mean sea level pressure

    */

    // loop over nodes
    for (int i=0; i<NNodes; i++) {
        
        double latitude = y[i]; 
        double longitude = x[i];

        // get values from the netcdf file
        double values[nFields];
        FluxesReader_global.GetScalars(longitude, latitude, values);
        // values contains the values above in the order we registered them
        // See above

        // rotate wind to cartesian grid
        double u_rot = Vx[i];
        double v_rot = Vy[i];
        if (rotate) {
          rotate_wind_fc(&longitude, &latitude, &Vx[i], &Vy[i], &Vz[i], &u_rot, &v_rot);
        }
        // DelU - the difference between wind and water currents
        delU_u[i] = values[0] - u_rot;
        delU_v[i] = values[1] - v_rot;

        // set up SST
        if (T[i] < 250.0) {
            // we're in degress C, so convert to K
            SST[i] = T[i] + 273.15;
        } else {
            // already in K - don't convert
            SST[i] = T[i];
        }
            
        // calculate parameters required for bulk forcing from ERA data
        z_data[i] = 2.0; // assume 2m height for t & q
        speed[i] = sqrt(pow(delU_u[i],2.0) + pow(delU_v[i],2.0));
        qs[i] = (1.0/air_density)*q1*exp(q2/SST[i]);
        avail[i] = true;
        t_2m[i] = values[7];
        double msl = values[8];
        // specific humidity from dewpoint and air temp
        // q = E_s/T*R
        // E_s vapor pressure in Pa
        // E_s = 611 * 10^(7.5Td/(273-Td)) where Td is dewpoint temp in deg C
        // T is air temp in K and R is gas constant for water vapor = 461.5 (J/kg*Kelvin)
        double temp = values[6] - 273.15;
        temp = 7.5*temp / (237.3+temp);
        temp = (610.78 * pow(10.0,temp));
        q[i] = 0.62197 * (temp / (msl + temp*(0.62197-1)));
        runoff[i] = values[4] / accumulated_correction;
        ppt[i] = values[5] / accumulated_correction;

        // fix integrated values - assume 6 hours
        //  - ssr
        solar[i] = (values[2] / accumulated_correction);
        //  - str
        thermal[i] = (values[3] / accumulated_correction);
    }
    delete [] x;
    delete [] y;
    delete [] z;
        
    // ncar_forcing requires that the inputs are array corresponding to
    // each point that a forcing is required. In order these are:
    //  - u_del - wind speed relative to currents (currents usually ignored)
    //  - t_2m - 2m temperature
    //  - ts - SST
    //  - q - specific humidity
    //  - qs - saturating humidity
    //  - z - height of point i
    //  - avail - array of booleans to say if data is available at point i
    //  - cd - these are the coefficients calculated by the sub routine
    //  - ch - that are then used to calculate QH, E and Tau from 
    //  - ce - Large and Yeager (2004), eq: 4a-d
    //  - ustar
    //  - bstar
    ncar_ocean_fluxes_fc(&NNodes,speed, t_2m, SST, q, qs, z_data, avail,
            cd, ch, ce, ustar, bstar);

   
    double heat_convert = ocean_density * ocean_heat_capacity;
    heat_convert = 1.0 / heat_convert;    
    double OneOverDensity = 1.0 / ocean_density;
    
    double E, Q_s, Q_l, Q_e, Q_h, Q_p, F;
    for (int i=0; i<NNodes; i++) {
        // from cd, ce and ch, calculate fluxes
        double tau_temp = OneOverDensity * air_density * cd[i] * speed[i];
        tau_u[i] = tau_temp * delU_u[i];
        tau_v[i] = tau_temp * delU_v[i];
        E = air_density * ce[i] * (q[i] - qs[i]) * speed[i]; // evap
        Q_s = solar[i] * (1.0-alpha);
        Q_l = thermal[i] - (sb * pow(SST[i],4.0)); //longwave
        Q_e = vapour_latent * E; //latent
        Q_h = air_density * air_specificHeat * ch[i] * (t_2m[i] - SST[i]) * speed[i]; //sensible
        Q_p = -fusion_latent * ppt[i];
        // E seems to be about a factor of a thousand out of the ERA40 data when using 
        // it for the freshwater fluxes
        // Dividing by ocean density appears to give the right answer...
        E = E / ocean_density;
        F = ppt[i] + E + runoff[i];

        Q_solar[i] = (solar[i]);
        Q_as[i] = heat_convert * (Q_s + Q_l + Q_e + Q_h + Q_p);
        F_as[i] = -1.0 * ref_salinity * F;
       
        
        // optional output
        if (e) {
          e[i] = E;
        }
        if (q_l) {
          q_l[i] = Q_l;
        }
        if (q_h) {
          q_h[i] = Q_h;
        }
        if (q_e) {
          q_e[i] = Q_e;
        }
        if (q_p) {
          q_p[i] = Q_p;
        }
    }

    //cout<<Q_as[0]<<endl;

  // clean up
  delete [] delU_u;
  delete [] delU_v;
  delete [] cd;
  delete [] ch;
  delete [] ce;
  delete [] t_2m;
  delete [] q;
  delete [] z_data;
  delete [] solar; 
  delete [] thermal;
  delete [] cloud;
  delete [] runoff;
  delete [] ppt;
  delete [] avail;
  delete [] SST;

}

#ifdef FORCINGERA40_UNITTEST
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

int main(int argc, char **argv){
    
    char *filename=argv[1];
    double accumulated_correction = 6.0*60.0*60.0;
    vtkXMLUnstructuredGridReader *reader = vtkXMLUnstructuredGridReader::New();
    reader->SetFileName(filename);
    reader->Update();
  
    vtkUnstructuredGrid *ug = reader->GetOutput();

    int NNodes = ug->GetNumberOfPoints();
    double time = 21601;  
    vector<double> T(NNodes, 0.0), S(NNodes, 35.0), Vx(NNodes, 0.0), Vy(NNodes, 0.0), Vz(NNodes, 0.0),
    F_as(NNodes, 0.0), Q_as(NNodes, 0.0), tau_u(NNodes, 0.0), tau_v(NNodes, 0.0),
    X(NNodes, 0.0), Y(NNodes, 0.0), Z(NNodes, 0.0), Q_solar(NNodes, 0.0), e(NNodes,0.0),
    q_l(NNodes,0.0), q_h(NNodes,0.0), q_e(NNodes, 0.0), q_p(NNodes, 0.0), x(NNodes,0.0), 
    y(NNodes,0.0), z(NNodes,0.0);
 
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
  
    // these are the data as simulated by ERA40
    FluxesReader FluxesReader_ERAdata;
    //FluxesReader_ERAdata.VerboseOn();
    FluxesReader_ERAdata.RegisterDataFile("1990_bats_in.nc");
    FluxesReader_ERAdata.SetSimulationTimeUnits("seconds since 1990-01-01 00:00:0");
    FluxesReader_ERAdata.SetTimeSeconds(time);
    FluxesReader_ERAdata.AddFieldOfInterest("sshf"); // 0 | surface sensible heat flux
    FluxesReader_ERAdata.AddFieldOfInterest("slhf"); // 1 | surface latent heat flux
    FluxesReader_ERAdata.AddFieldOfInterest("ssr");  // 2 | surface solar radiation
    FluxesReader_ERAdata.AddFieldOfInterest("str");  // 3 | surface thermal radiation
    FluxesReader_ERAdata.AddFieldOfInterest("ewss"); // 4 | east-west surface stress 
    FluxesReader_ERAdata.AddFieldOfInterest("nsss"); // 5 | north-south surface stress 
    FluxesReader_ERAdata.AddFieldOfInterest("e");    // 6 | evaporation

    // these data are the input to our fluxes routine
    FluxesReader_global.RegisterDataFile("1990_bats_out.nc");
    //FluxesReader_global.VerboseOn();
    FluxesReader_global.SetSimulationTimeUnits("seconds since 1990-01-01 00:00:0.0");
    FluxesReader_global.SetTimeSeconds(time);
    FluxesReader_global.AddFieldOfInterest("10u");  //  0   | 10 metre U wind component
    FluxesReader_global.AddFieldOfInterest("10v");  //  1   | 10 metre V wind component
    FluxesReader_global.AddFieldOfInterest("ssrd"); //  2   | Surface solar radiation downwards
    FluxesReader_global.AddFieldOfInterest("strd"); //  3   | Surface thermal radiation downwards
    FluxesReader_global.AddFieldOfInterest("ro");   //  4   | Runoff
    FluxesReader_global.AddFieldOfInterest("tp");   //  5   | Total precipitation
    FluxesReader_global.AddFieldOfInterest("2d");   //  6   | Dewpoint temp at 2m 
    FluxesReader_global.AddFieldOfInterest("2t");   //  7   | Air temp at 2m 
    FluxesReader_global.AddFieldOfInterest("msl");  //  8   | Mean sea level pressure 

    int ret;
    int const nFields_in=9, nFields_out=7;
    double values_in[nFields_in];
    double values_out[nFields_out];

/**
    // get the global SST from the "output" ERA40 data
    ret = projections(NNodes, &x[0], &y[0], &z[0], "cart", "spherical");
    if (ret != 0) {
        cerr<<"Error converting coord system"<<endl;
    }
    
    // calculate SST for each node
    for (int i=0; i<NNodes; i++) {
      FluxesReader_global.GetScalars(x[i], y[i], values_in);
      FluxesReader_ERAdata.GetScalars(x[i], y[i], values_out);
      double temp =  (values_in[3] - values_out[3]);
      temp = temp / (6.0 * 60.0 * 60.0); // correcting for accumulated values
      temp = temp / 5.67e-8;             // SB constant
      temp = pow(temp, 1.0/4.0);         // finally, get temperature
      T[i] = temp;
    }

    // now calculate fluxes for the whole globe and put on the output VTU file for
    // visualisation
    get_era40_fluxes_fc(&time, &X[0], &Y[0], &Z[0], &T[0], &Vx[0], &Vy[0], &Vz[0], &F_as[0],
                        &Q_as[0], &tau_u[0], &tau_v[0], &Q_solar[0], &NNodes, true, &q_l[0], &q_h[0],
                        &e[0]);

    // fix the values so they match ERA40 values straight from the NetCDF files.
    for (int i=0; i<NNodes; i++) {
      q_l[i] *= accumulated_correction;
      q_h[i] *= accumulated_correction;
      e[i] *= accumulated_correction;
      tau_u[i] *= accumulated_correction;
      tau_v[i] *= accumulated_correction;
    }

    // add to VTU file
    vtk_add_scalar(F_as, "F_as", ug);
    vtk_add_scalar(Q_as, "Q_as", ug);
    vtk_add_scalar(Q_solar, "Solar Radiation", ug);
    vtk_add_scalar(tau_u, "tau_u", ug);
    vtk_add_scalar(tau_v, "tau_v", ug);
    vtk_add_scalar(T, "SST", ug);
    vtk_add_scalar(q_l, "Latent Heat Flux", ug);
    vtk_add_scalar(q_h, "Sensible Heat Flux", ug);
    vtk_add_scalar(e, "Evaporation", ug);
    vtkXMLUnstructuredGridWriter *writer = vtkXMLUnstructuredGridWriter::New();
    writer->SetFileName("fluxes.vtu");
    writer->SetInput(ug);
    writer->Write();
   **/ 

    // now do a single node so we can compare values
    // grab the data and construct the input variables we would otherwise get
    // from ICOM. We overwrite the global data at point [0] here
    // convert from Cart to long-lat
    int n = 1;
    // Ocean Station Papa
    // -3.35835e+06,-2.35154e+06,4.88594e+06
    //X[0] = -3.35835e+06;
    //Y[0] = -2.35154e+06;
    //Z[0] = 4.88594e+06;
    //Bermuda
    X[0]=2775060;
    Y[0]=-5731974;
    Z[0]=352332;
    y[0] = Y[0] , x[0] = X[0] , z[0] = Z[0];
    ret = projections(n, &x[0], &y[0], &z[0], "cart", "spherical");
    if (ret != 0) {
        cerr<<"Error converting coord system"<<endl;
    }
    
    double Tau_u_out, Tau_v_out, F_out, Q_out;
    FluxesReader_ERAdata.GetScalars(x[0], y[0], values_out);
    FluxesReader_global.GetScalars(x[0], y[0], values_in);
    double temp =  (values_in[3] - values_out[3]);
    temp = temp / accumulated_correction; // correcting for accumulated values
    temp = temp / 5.67e-8;                // SB constant
    temp = pow(temp, 1.0/4.0);
    T[0] = temp;
    double ocean_density = 1027.0;
    double ref_salinity = 35.0;

    // now let's work out the fluxes for our test location
    Tau_u_out = values_out[4]/(accumulated_correction*ocean_density);
    Tau_v_out = values_out[5]/(accumulated_correction*ocean_density);
    F_out = values_in[5] + values_out[6] + values_in[4];
    F_out = F_out / (accumulated_correction * ocean_density) * ref_salinity;
    double Q_p = -3.337e5 * values_in[5];
    Q_out = values_out[2] + values_out[3] + values_out[0] + values_out[1] + Q_p;
    // 4000 = ocean heat capacity
    Q_out = Q_out / (accumulated_correction*ocean_density*4000.0);

    // get fluxes for test location
    get_era40_fluxes_fc(&time, &X[0], &Y[0], &Z[0], &T[0], &S[0], &Vx[0], &Vy[0], &Vz[0], &F_as[0],
                        &Q_as[0], &tau_u[0], &tau_v[0], &Q_solar[0], &n, false, &q_l[0], &q_h[0], &e[0],
                        &q_e[0], &q_p[0]);

    cout <<endl<<"Fluxes for "<<x[0]<<","<<y[0]<<endl<<endl;

    cout <<"Flux                 |    ERA40 value     |  Our value   "<<endl;
    cout <<"---------------------+--------------------+---------------"<<endl;
    printf("U stress (m2/s2)     | %-18g | %-18g\n",Tau_u_out,tau_u[0]);
    printf("V stress (m2/s2)     | %-18g | %-18g\n",Tau_v_out,tau_v[0]);
    printf("heat flux (K/s)      | %-18g | %-18g\n",Q_out,Q_as[0]);
    printf("Salinity flux        | %-18g | %-18g\n",F_out,F_as[0]);
    printf("Sens heat flux (W/m2)| %-18g | %-18g\n",values_out[0]/accumulated_correction,q_h[0]);
    printf("Lat heat flux (W/m2) | %-18g | %-18g\n",values_out[1]/accumulated_correction,q_e[0]);
    printf("Evaporation (m/s)    | %-18g | %-18g\n",values_out[6]/accumulated_correction,e[0]);
    printf("Thermal (W/m2)       | %-18g | %-18g\n",values_out[3]/accumulated_correction,q_l[0]);
    printf("PPT heat (W/m2)      | %-18g | %-18g\n",Q_p/accumulated_correction,q_p[0]);
    printf("Solar rad. (W/m2)    | %-18g | %-18g\n",values_out[2]/accumulated_correction,values_in[2]/accumulated_correction*(1-0.066));
    cout<<endl<<endl;

    cout <<" Time (s),ERA40,Our Q,PPT,E,ICOM E,e-p ICOM,e-p ERA40"<<endl;
    // let's print out heat flux over 28 days, every 6 hours to compare
    for (int t = 21600 ; t < 31536000; t += 86400) {
        FluxesReader_ERAdata.SetTimeSeconds(t);
        FluxesReader_ERAdata.GetScalars(x[0], y[0], values_out);
        FluxesReader_global.SetTimeSeconds(t);
        FluxesReader_global.GetScalars(x[0], y[0], values_in);
        Tau_u_out = values_out[4]/(accumulated_correction*ocean_density);
        Tau_v_out = values_out[5]/(accumulated_correction*ocean_density);


        double Q_p = -3.337e5 * values_in[5];
        Q_out = values_out[2] + values_out[3] + values_out[0] + values_out[1] + Q_p;
        // 4000 = ocean heat capacity
        Q_out = Q_out / (accumulated_correction*ocean_density*4000.0);
        double temp =  (values_in[3] - values_out[3]);
        temp = temp / accumulated_correction; // correcting for accumulated values
        temp = temp / 5.67e-8;              // SB constant
        temp = pow(temp, 1.0/4.0);
        T[0] = temp;


        double t2 = t;
        get_era40_fluxes_fc(&t2, &X[0], &Y[0], &Z[0], &T[0], &S[0], &Vx[0], &Vy[0], &Vz[0], &F_as[0],
                        &Q_as[0], &tau_u[0], &tau_v[0], &Q_solar[0], &n, true, &q_l[0], &q_h[0], &e[0],
                        &q_e[0], &q_p[0]);

        printf("%d,%g,%g,%g,%g,%g,%g,%g\n",t,Q_out,Q_as[0],values_in[5]/accumulated_correction,values_out[6]/accumulated_correction,e[0],F_as[0],(values_in[5]/accumulated_correction+values_out[6]/accumulated_correction)*-35.);

    }

    
}
#endif

#ifdef WIND
#include <iostream>
#include <fstream>
#include <vector>

int main(int argc, char **argv){

    int NNodes = 1;

    int start = 0;
    int end = 3*365*24*60*60;
    int step = 60*60;
    
    // assume temperature is an even 280K
    vector<double> T(NNodes, 280.0), Vx(NNodes, 0.0), Vy(NNodes, 0.0), Vz(NNodes, 0.0),
    F_as(NNodes, 0.0), Q_as(NNodes, 0.0), tau_u(NNodes, 0.0), tau_v(NNodes, 0.0),
    X(NNodes, 0.0), Y(NNodes, 0.0), Z(NNodes, 0.0), Q_solar(NNodes, 0.0), e(NNodes,0.0),
    q_l(NNodes,0.0), q_h(NNodes,0.0), q_e(NNodes, 0.0), q_p(NNodes, 0.0), x(NNodes,0.0), 
    y(NNodes,0.0), z(NNodes,0.0);

    // Ocean Station Papa
    X[0] = -3070980.0;
    Y[0] = -2150320.0;
    Z[0] = 5160020.0;
    y[0] = Y[0] , x[0] = X[0] , z[0] = Z[0];
    int ret = projections(NNodes, &x[0], &y[0], &z[0], "cart", "spherical");
    if (ret != 0) {
        cerr<<"Error converting coord system"<<endl;
    }

    // these data are the input to our fluxes routine
    FluxesReader_global.RegisterDataFile("1990_bats.nc");
    //FluxesReader_global.VerboseOn();
    FluxesReader_global.SetSimulationTimeUnits("seconds since 1990-01-01 00:00:0.0");
    FluxesReader_global.AddFieldOfInterest("tcc");  //  0   | Total Cloud Cover
    FluxesReader_global.AddFieldOfInterest("10u");  //  1   | 10 metre U wind component
    FluxesReader_global.AddFieldOfInterest("10v");  //  2   | 10 metre V wind component
    FluxesReader_global.AddFieldOfInterest("ssrd"); //  3   | Surface solar radiation downwards
    FluxesReader_global.AddFieldOfInterest("strd"); //  4   | Surface thermal radiation downwards
    FluxesReader_global.AddFieldOfInterest("ro");   //  5   | Runoff
    FluxesReader_global.AddFieldOfInterest("tp");   //  6   | Total precipitation
    FluxesReader_global.AddFieldOfInterest("q");    //  7   | Specific humidity at z
    FluxesReader_global.AddFieldOfInterest("2t");   //  8   | Air temp at 2m 

    // open file to dump to
    ofstream output;
    output.open ("wind_stress.csv");

    // loop over time, grabbing fluxes, printing wind stress, every hour.
    for( double time=start; time<end; time=time+step) {
        get_era40_fluxes_fc(&time, &X[0], &Y[0], &Z[0], &T[0], &S[0], &Vx[0], &Vy[0], &Vz[0], &F_as[0],
                        &Q_as[0], &tau_u[0], &tau_v[0], &Q_solar[0], &NNodes, true, &q_l[0], &q_h[0],
                        &e[0]);

        // write time, u, v
        output<< time <<","<<tau_u[0]<<","<<tau_v[0]<<endl;
    }

    // close file
    output.close();
}

#endif
