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
#include "Particle.h"

Particle::Particle(double coordinates[], int dim, int cell, double local_coords[]) :
  position(dim), local_coords(dim+1)
{
  this->cell = cell;

  /* Deep copy position coordinates */
  for (int i=0; i<dim; i++) this->position[i] = coordinates[i];
  for (int i=0; i<dim+1; i++) this->local_coords[i] = local_coords[i];
}

Particle::~Particle() {}

void Particle::view()
{
  cout << "Static Particle :: " << endl;
  cout << "\tPosition: " << position[0] << ", " << position[1];
  cout << "  \tcell: " << cell << endl << "\tLocal coords: ";
  cout << local_coords[0] << ", " << local_coords[1] << ", " << local_coords[2];
  cout << endl;
}
