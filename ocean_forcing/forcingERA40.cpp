/*  Copyright (C) 2009 Imperial College London and others.

    Please see the AUTHORS file in the main source directory for a full list
    of copyright holders.

    Prof. C Pain
    Applied Modelling and Computation Group
    Department of Earth Science and Engineering
    Imperial College London

    amcgsoftware@imperial.ac.uk

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
#include "BulkForcing.h"
#include "FluxesReader.h"
#include "spud"
#include "global_parameters.h"
#include "coordinates.h"

using namespace std;
using namespace Spud;

void get_era40_fluxes_fc(double *time, const double *X, const double *Y, const double *Z,
                     double *T, const double *Vx, const double *Vy, const double *Vz, double *salinity,
                     double *F, double *Q, double *tau_u, double *tau_v, double *solar,
                     const int *n, bool rotate, int *bulk_formula ) {

    // physical constants
    const double air_density = 1.22,
               q1 = 0.98*640380,
               q2 = -5107.4,
               kelvin_centrigrade = 273.15;

    // Set up for ERA40 data. Options to change rate are set below
    double accumulated_correction = 6.0*60.0*60.0; // Assumes data every 6 hours.

    // problem constants
    const int nFields = 9;
    const int NNodes = *n;

    double *delU_u = new double[NNodes];
    double *delU_v = new double[NNodes];
    double *t_2m = new double[NNodes];
    double *q = new double[NNodes];
    double *z_data = new double[NNodes];
    double *Q_solar = new double[NNodes];
    double *thermal = new double[NNodes];
    double *cloud = new double[NNodes];
    double *runoff = new double[NNodes];
    double *ppt = new double[NNodes];
    double *speed = new double[NNodes];
    double *qs = new double[NNodes];
    double *SST = new double[NNodes];

    FluxesReader_global.SetTimeSeconds(*time);

    if(have_option("/ocean_forcing/bulk_formulae/input_file_type/type::ERA40/no_accumulation")) {
        cout<<"Found flag"<<endl;
        accumulated_correction = 1;
    }


    double surface_radius = get_surface_radius();
    double longitude[NNodes];
    double latitude[NNodes];
    double height[NNodes];
    double u_rot = 0.0;
    double v_rot = 0.0;
    double w_rot = 0.0;

    // loop over nodes
    for (int i=0; i<NNodes; i++) {

        // Rotate ocean surface velocity to zonal-meridional-vertical. Also
        //  Transforms cartesian position components into lon-lat-height.
        double u_cart = Vx[i];
        double v_cart = Vy[i];
        double w_cart = Vz[i];
        double x_cart = X[i];
        double y_cart = Y[i];
        double z_cart = Z[i];
        if (rotate) {
          vector_cartesian_2_lon_lat_height_c(&u_cart, &v_cart, &w_cart,
                                              &x_cart, &y_cart, &z_cart,
                                              &u_rot, &v_rot, &w_rot,
                                              &longitude[i], &latitude[i], &height[i],
                                              &surface_radius);
        }else{
          cartesian_2_lon_lat_height_c(&x_cart, &y_cart, &z_cart,
                                       &longitude[i], &latitude[i], &height[i],
                                       &surface_radius);
        }

        /*
         *Get values from the netcdf file.
         * The "values" array below contains the scalar
         * values from the ERA40 netCDF file. The table
         * below shows which ERA40 parameter is at which
         * array index, along with a physical meaning
         *
        ERA40 field | Index | Physical meaning
        ------------+-------+-----------------
        u10         |   0   | 10 metre U wind component
        v10         |   1   | 10 metre V wind component
        ssrd        |   2   | Surface solar radiation
        strd        |   3   | Surface thermal radiation
        ro          |   4   | Runoff
        tp          |   5   | Total precipitation
        d2m         |   6   | Dewpoint temp at 2m
        t2m         |   7   | Air temp at 2m
        msl         |   8   | Mean sea level pressure
        */
        double values[nFields];
        FluxesReader_global.GetScalars(longitude[i], latitude[i], values);

        // DelU - the difference between wind and water currents
        delU_u[i] = values[0] - u_rot;
        delU_v[i] = values[1] - v_rot;

        // set up SST
        if (T[i] < 250.0) {
            // we're in degress C, so convert to K
            SST[i] = T[i] + kelvin_centrigrade;
        } else {
            // already in K - don't convert
            SST[i] = T[i];
        }

        // calculate parameters required for bulk forcing from ERA data
        z_data[i] = 2.0; // assume 2m height for t & q
        speed[i] = sqrt(pow(delU_u[i],2.0) + pow(delU_v[i],2.0));
        qs[i] = (1.0/air_density)*q1*exp(q2/SST[i]); // Sea Surface specific Humidity
        t_2m[i] = values[7];
        double msl = values[8];
        // specific humidity (q) from air temp
        double temp = values[6] - kelvin_centrigrade;
        temp = 7.5*temp / (237.3+temp);
        temp = (610.78 * pow(10.0,temp));
        q[i] = 0.62197 * (temp / (msl + temp*(0.62197-1)));
        // Alternative formula for q from Fairall (COARE 3)
        // Appears to give better Tau, but much worse Q and F, when using NCAR param.
        //temp = 6.112*exp(17.502*temp/(temp+241.0))*(1.0007+3.46e-6*msl);
        //q[i] = temp*622./(msl-.378*temp);
        runoff[i] = values[4] / accumulated_correction;
        ppt[i] = values[5] / accumulated_correction;

        // fix integrated values - assume 6 hours
        //  - ssr
        solar[i] = (values[2] / accumulated_correction);

        //  - str
        thermal[i] = (values[3] / accumulated_correction);
    }


    switch (*bulk_formula) {

        case COARE3:
            coare_forcing_c_fc(&NNodes, speed, t_2m, SST, q, qs, delU_u, delU_v,
                    ppt, runoff, salinity, thermal, solar, Q_solar, Q, F, tau_u, tau_v);
            break;
        case KARA:
            kara_forcing_c_fc(&NNodes, speed, t_2m, SST, q, qs, delU_u, delU_v,
                    ppt, runoff, salinity, thermal, solar, Q_solar, Q, F, tau_u, tau_v);
            break;

        case NCAR:
        default  :
            ncar_forcing_c_fc(&NNodes, speed, t_2m, SST, q, qs, delU_u, delU_v,
                    ppt, runoff, salinity, thermal,  solar, Q_solar, Q, F, tau_u, tau_v);
    }

    // clean up
    delete [] delU_u;
    delete [] delU_v;
    delete [] t_2m;
    delete [] q;
    delete [] z_data;
    delete [] Q_solar;
    delete [] thermal;
    delete [] cloud;
    delete [] runoff;
    delete [] ppt;
    delete [] SST;
    delete [] speed;
    delete [] qs;

}
