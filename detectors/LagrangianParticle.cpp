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
#include "LagrangianParticle.h"

extern "C" {
  /* The following routines should eventually become dynamically set callbacks */
  void get_local_coordinates(int dim, int element, double *coordinates,
                             double *local_coordinates);

  int get_element_neighbour(int element, int neigh_i);
}

LagrangianParticle::LagrangianParticle(double coordinates[], int dim, int cell,
                                       double local_coords[], char *name) :
  Particle(coordinates, dim, cell, local_coords, name), update_vector(dim)
{
  for (size_t i=0; i<this->position.size(); i++) {
    this->update_vector[i] = this->position[i];
  }
}

LagrangianParticle::~LagrangianParticle() {}

void LagrangianParticle::view()
{
  cout << "Lagrangian Particle :: " << endl;
  cout << "\tPosition: " << position[0] << ", " << position[1];
  cout << "  \tcell: " << cell << endl << "\tLocal coords: ";
  cout << local_coords[0] << ", " << local_coords[1] << ", " << local_coords[2];
  cout << endl;
}

int    LagrangianParticle::rk4_stages = 4;
double LagrangianParticle::rk4_timestep_weights[4] = { 1./6., 1./3., 1./3., 1./6. };
double LagrangianParticle::rk4_stage_matrix[4][4] = { {0., 0., 0., 0.},
                                                      {.5, 0., 0., 0.},
                                                      {0., .5, 0., 0.},
                                                      {0., 0., 1., 0.} };

void LagrangianParticle::rk4_advection(void *velocity_field, double dt)
{
  double **k = new double*[rk4_stages];
  double *update_lcoords = new double[local_coords.size()];

  for (int stage=0; stage<rk4_stages; stage++) {
    k[stage] = new double[position.size()];

    /* Stage vector is computed by evaluating velocity at the current position */
    get_local_coordinates(update_vector.size(), cell,
                          update_vector.data(), update_lcoords);
    vector_field_value(velocity_field, update_lcoords, k[stage]);

    if (stage < rk4_stages-1) {
      /* Update vector maps from current position to place required
         for computing next stage vector */
      for (size_t i=0; i<this->position.size(); i++) update_vector[i] = position[i];
      for (int s=0; s<=stage; s++) {
        for (size_t i=0; i<this->update_vector.size(); i++) {
          update_vector[i] += dt * rk4_stage_matrix[stage+1][s] * k[s][i];
        }
      }
    } else {
      /* Update vector maps from current position to final position */
      for (size_t i=0; i<this->position.size(); i++) update_vector[i] = position[i];
      for (int s=0; s<rk4_stages; s++) {
        for (size_t i=0; i<this->update_vector.size(); i++) {
          update_vector[i] += dt * rk4_timestep_weights[s] * k[s][i];
        }
      }
    }

    update_parametic(update_vector, 1.e-10);
  }
  for (size_t i=0; i<this->position.size(); i++) position[i] = update_vector[i];

  for (int stage=0; stage<rk4_stages; stage++) delete [] k[stage];
  delete [] k;
  delete [] update_lcoords;
}

void LagrangianParticle::update_parametic(vector<double> coordinates, double tolerance)
{
  int neigh_i, dim = this->local_coords.size();
  double *lcoords = new double[dim];
  double *min_lcoord;

  while (*min_lcoord < -tolerance) {
    get_local_coordinates(coordinates.size(), cell, coordinates.data(), lcoords);
    min_lcoord = min_element(lcoords, lcoords+dim);

    neigh_i = min_lcoord - lcoords;
    this->cell = get_element_neighbour(this->cell, neigh_i);
  }
  delete [] lcoords;
}
