/*  Copyright (C) 2006 Imperial College London and others.
    
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
#ifndef BULKFORCING_H
#define BULKFORCING_H

#include <string>
#include <stdio.h>
#include <stdlib.h>


#define NCAR 0
#define COARE3 1
#define KARA 2


extern "C" {
#define ncar_forcing_c_fc F77_FUNC(ncar_forcing_c, NCAR_FORCING_C)
  extern void ncar_forcing_c_fc(const int *points, double *speed, double *air_temp, 
                              double *sst, double *spec_humidity, double *sea_surface_humidity, 
                              double *U, double *V, double *ppt, double *runoff, double *salinity,
                              double *solar, double *thermal, double *Q_solar, double *Q, double *F, 
                              double *tau_u, double *tau_v);

#define coare_forcing_c_fc F77_FUNC(coare_forcing_c, COARE_FORCING_C)
  extern void coare_forcing_c_fc(const int *points, double *speed, double *air_temp, 
                              double *sst, double *spec_humidity, double *sea_surface_humidity, 
                              double *U, double *V, double *ppt, double *runoff, double *salinity,
                              double *solar, double *thermal, double *Q_solar, double *Q, double *F, 
                              double *tau_u, double *tau_v);

#define kara_forcing_c_fc F77_FUNC(kara_forcing_c, KARA_FORCING_C)
  extern void kara_forcing_c_fc(const int *points, double *speed, double *air_temp, 
                              double *sst, double *spec_humidity, double *sea_surface_humidity, 
                              double *U, double *V, double *ppt, double *runoff, double *salinity,
                              double *solar, double *thermal, double *Q_solar, double *Q, double *F, 
                              double *tau_u, double *tau_v);

#define rotate_wind_fc F77_FUNC(rotate_wind, ROTATE_WIND)
  extern void rotate_wind_fc(double *longitude, double *latitude, const double *r3u, 
                            const double *r3v, const double *r3w,
                            double *u, double *v);


#define get_era40_fluxes_fc F77_FUNC(get_era40_fluxes, GET_ERA40_FLUXES)
    void get_era40_fluxes_fc(double *time, const double *X, const double *Y, const double *Z, 
                     double *T, const double *Vx, const double *Vy, const double *Vz, double *Sal,
                     double *F_as, double *Q_as, double *tau_u, double *tau_v, double *Q_solar,
                     const int *NNodes, bool rotate, int *bulk_formula);
}


#endif
