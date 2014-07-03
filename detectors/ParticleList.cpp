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
#include "ParticleList.h"

extern "C" {
  void find_enclosing_cell(double *coordinates, int dim, int *cell, double *lcoords);
}

ParticleList::ParticleList() {}

ParticleList::~ParticleList()
{
  plist.clear();
}

void ParticleList::view()
{
  cout << "ParticleList :: size " << plist.size() << endl;
  for (list<Particle*>::iterator it = plist.begin(); it != plist.end(); ++it) {
    (*it)->view();
  }
}

void ParticleList::add_particle(Particle *p) {
  plist.push_back(p);
}

list<Particle*>::iterator ParticleList::begin()
{
  return plist.begin();
}

list<Particle*>::iterator ParticleList::end()
{
  return plist.end();
}

void ParticleList::advect_lagrangian(void *velocity_field, double dt) {
  for (list<Particle*>::iterator p = plist.begin(); p != plist.end(); ++p) {
    LagrangianParticle *lagp = dynamic_cast<LagrangianParticle*>(*p);
    if (lagp != NULL) {
      lagp->rk4_advection(velocity_field, dt);
    }
  }
}

extern "C" {
  ParticleList *create_particle_list()
  {
    return new ParticleList();
  }

  void particle_list_view(ParticleList *plist)
  {
    plist->view();
  }

  void particle_list_new_particle(ParticleList *plist, double coordinates[],
                                  int dim, char *name)
  {
    /* Establish cell and local coordinates via callback */
    double *lcoords = new double[dim+1];
    int cell = 0;
    find_enclosing_cell(coordinates, dim, &cell, lcoords);

    if (cell > 0) {
      Particle *p = new Particle(coordinates, dim, cell, lcoords, name);
      plist->add_particle(p);
    }
    delete [] lcoords;
  }

  void particle_list_new_lagrangian_particle(ParticleList *plist, double coordinates[],
                                             int dim, char *name)
  {
    /* Establish cell and local coordinates via callback */
    double *lcoords = new double[dim+1];
    int cell = 0;
    find_enclosing_cell(coordinates, dim, &cell, lcoords);

    if (cell > 0) {
      LagrangianParticle *p = new LagrangianParticle(coordinates, dim, cell,
                                                     lcoords, name);
      plist->add_particle(p);
    }
    delete [] lcoords;
  }

  void particle_list_advect_lagrangian(ParticleList *plist, void *velocity_field,
                                       double dt) {
    plist->advect_lagrangian(velocity_field, dt);
  }
}
