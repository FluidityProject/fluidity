/* Copyright (C) 2006 Imperial College London and others.

Please see the AUTHORS file in the main source directory for a full list
of copyright holders.

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
//
// Explination of sam options
//

// options[0]:  m // Target number of partitions. If this is 0, then
//                // the number of ranks in MPI_COMM_WORLD is assumed.
// options[1]:  1 // clean repartitioning into nparts
//              2 // re-partitioning (diffusive) - local optimization
//              3 // re-partitioning (directed diffusion) - global optimization
//              4 // re-partitioning (repartition-remapping)
// options[2]:  1 // homogenious processing power
//              2 // hetrogenious processing power
// options[3]:  1 // no node weights
//              2 // calculate edge weights based on the projected 
//                   density of nodes in the future
// options[4]:  1 // no edge weights
//              2 // calculate edge weights optimizing for adaptivity
//              3 // calculate edge weights based on curvature of fields
// options[5]:  1 // simple mesh
//              2 // mixed formulation
#include "confdefs.h"

#include <map>
#include <string>
#include <vector>

#include "c++debug.h"
#include "Mesh.h"
#include "samtypes.h"

using namespace std;

#define sam_fc F77_FUNC(sam, SAM)
#define sam_init_fc F77_FUNC(sam_init_c, SAM_INIT_C)
#define sam_migrate_fc F77_FUNC(sam_migrate_c, SAM_MIGRATE_C)
#define sam_query_fc F77_FUNC(sam_query_c, SAM_QUERY_C)
#define sam_export_mesh_fc F77_FUNC(sam_export_mesh_c, SAM_EXPORT_MESH_C)
#define sam_export_halo_fc F77_FUNC(sam_export_halo_c, SAM_EXPORT_HALO_C)
#define sam_export_phalo_fc F77_FUNC(sam_export_phalo_c, SAM_EXPORT_PHALO_C)
#define sam_cleanup_fc F77_FUNC(sam_cleanup_c, SAM_CLEANUP_C)
#define sam_add_field_fc F77_FUNC(sam_add_field_c, SAM_ADD_FIELD_C)
#define sam_pop_field_fc F77_FUNC(sam_pop_field_c, SAM_POP_FIELD_C)
#define sam_export_node_ownership_fc F77_FUNC(sam_export_node_ownership_c, SAM_EXPORT_NODE_OWNERSHIP_C)

// "Mesh" is the principle object that does everything for SAM.
Mesh *mesh = NULL;
vector<int> decomp_opts;
vector<int> noddom;
  
extern "C"{  
  void sam_init_fc(
                   const int* dim, const int *NNodes,  const int *NElems, const int *NSElems, 
                   const int GATHER[], const int ATOSEN[],
                   const int SCATER[], const int ATOREC[],
                   const int* ncolga, const int* nscate, const int* nprocs,
                   const int ENLIST[], const int *nloc,
                   const int SNLIST[], const int SURFID[], const int *snloc,
                   const samfloat_t NODX[], const samfloat_t NODY[], const samfloat_t NODZ[],
                   samfloat_t Metric[],     const samfloat_t FIELDS[], const int *NFIELDS, 
                   const int options[], const samfloat_t *fxnl_tol){
#ifdef HAVE_MPI
    if(mesh != NULL)
    {
      delete mesh;
      mesh = NULL;
    }
    mesh = new Mesh;

    CHECK( options[5] );
    mesh->mixed_formulation(options[5] == 2);
  
  
    if( mesh->mixed_formulation() )
      ECHO("Mixed-formulation capabilities have been enabled.");
  
    // Slurp in Fluiditys' mesh
    ECHO("Importing mesh");
    mesh->import_fluidity(*dim, *NNodes, *NElems, *NSElems, *NNodes,
			 *NFIELDS, NODX, NODY, NODZ, Metric, FIELDS,
			 ENLIST, *nloc,
			 SNLIST, SURFID, *snloc,
			 ATOREC, SCATER,
			 ATOSEN, GATHER);
    ECHO("Mesh imported");

    // Build options.
    decomp_opts.resize(5);
    for(int i=0;i<5;i++)
      decomp_opts[i] = options[i];

    mesh->set_functional_tol(*fxnl_tol);
  
    if( (options[2] == 2)&&(options[3] == 2) ){
      cout<<"WARNING: It has been requested to load balance the mesh on a system with "          
          << "heterogeneous point-to-point bandwidth and heterogeneous processing power."
          << "This is a non-trivial request. SAM promises to do his best, but don't "     
          << "expect wonders."<< endl;
    }
#endif
    return;
  }

  void sam_migrate_fc(void)
  {
#ifdef HAVE_MPI
    assert(mesh != NULL);
    
    // Graph partitioning.
    ECHO("Partitioning mesh...");
    noddom = mesh->decomp( decomp_opts );
    ECHO("...partitioned.");
  
    // Migrate Mesh according to noddom.
    ECHO("Migrating mesh...");
    mesh->migrate( noddom );
    ECHO("...migrated.");
  
    if(mesh->mixed_formulation()){
      mesh->invent_pressure_mesh();
      mesh->formHalo2();
    }
#endif
    return;
  }

  void sam_query_fc(int* NONODS, int* TOTELE, int* STOTEL,
                    int* ncolga, int* nscate, int* pncolga, int* pnscate)
  {
#ifdef HAVE_MPI
    assert(mesh != NULL);

    mesh->do_element_headcount();
    *NONODS = mesh->node_list.size();
    *TOTELE = mesh->num_elements("volume");
    *STOTEL = mesh->num_elements("surface");

    *ncolga = mesh->get_ncolga();
    *nscate = mesh->get_nscate();

    if(mesh->mixed_formulation())
    {
      *pncolga = mesh->get_pncolga();
      *pnscate = mesh->get_pnscate();
    }
    else
    {
      *pncolga = -1;
      *pnscate = -1;
    }
#endif
    return;
  }

  void sam_export_mesh_fc(int* nonods, int* totele, int* stotel, int* nloc, int* snloc,
      samfloat_t* NODX, samfloat_t* NODY, samfloat_t* NODZ, int* ENLIST, int* SENLIST, int* SURFID)
  {
#ifdef HAVE_MPI
    assert(mesh != NULL);

    // Positions field
    {
      int i = 0;
      for(deque<Node>::iterator in=mesh->node_list.begin(); in != mesh->node_list.end(); ++in){
        switch(in->get_size_x()){
	  case 3:
            NODZ[i] = in->get_z();
	  case 2:
            NODY[i] = in->get_y();
	  case 1:
            NODX[i] = in->get_x();
            break;
	  default:
	    ERROR("Invalid dimension");
	}
        i++;
      }
    }
    
    int num_elements = mesh->element_list.size();
    vector<int> elm_owner(num_elements);
    for(int i=0;i<num_elements;i++){
      vector<unn_t> enl(mesh->element_list[i].get_enlist());
      
      vector<unn_t>::iterator it=enl.begin();
      int gnn = mesh->unn2gnn(*it);
      elm_owner[i] = mesh->node_list[gnn].get_future_owner();
      
      for(;it!=enl.end();++it){
        gnn = mesh->unn2gnn(*it);
        elm_owner[i] = min(elm_owner[i], (int)mesh->node_list[gnn].get_future_owner());
      }
    }

    int MyRank = MPI::COMM_WORLD.Get_rank();

    { // Compress volume element-node lists
      int i = 0;
      vector<int> halo_elements;
      for(int e=0;e<num_elements;e++){
        unsigned char type = mesh->element_list[e].get_flags();
        if( type & ELM_VOLUME ){
          vector<unn_t> enl(mesh->element_list[e].get_enlist());
                    
          if(elm_owner[e]==MyRank){
            for(vector<unn_t>::const_iterator it=enl.begin(); it!=enl.end(); ++it){
              ENLIST[i++] = mesh->unn2gnn(*it) + 1;
            }
          }else{
            for(vector<unn_t>::const_iterator it=enl.begin(); it!=enl.end(); ++it){
              halo_elements.push_back(mesh->unn2gnn(*it) + 1);
            }
          }
        }
      }
      for(vector<int>::const_iterator it=halo_elements.begin();it!=halo_elements.end();++it){
        ENLIST[i++]=*it;
      }
    }

    { 
      // Compress surface element-node list and write the surface id's.    
      mesh->do_element_headcount();
#ifndef NDEBUG
      int NewNSElems = mesh->num_elements("surface");
#endif
      int pos=0;
      int i = 0;
      vector<int> halo_elements, halo_element_ids;
      for(int e=0;e<num_elements;e++){
        unsigned char type = mesh->element_list[e].get_flags();
        if( type & ELM_SURFACE ){
          assert(pos<NewNSElems);
          
          vector<unn_t> enl( mesh->element_list[e].get_enlist() );
          const vector<int>& ifields = mesh->element_list[e].get_ifields();
          if(elm_owner[e]==MyRank){
            for(vector<unn_t>::const_iterator it = enl.begin(); it != enl.end(); ++it){
              SENLIST[i++] = mesh->unn2gnn(*it) + 1;
            }
            SURFID[pos++] = ifields[0];
          }else{
            for(vector<unn_t>::const_iterator it = enl.begin(); it != enl.end(); ++it){
              halo_elements.push_back(mesh->unn2gnn(*it) + 1);
            }
            halo_element_ids.push_back(ifields[0]);            
          }
        }
      }
      for(vector<int>::const_iterator it=halo_elements.begin();it!=halo_elements.end();++it){
        SENLIST[i++]=*it;
      }
      for(vector<int>::const_iterator it=halo_element_ids.begin();it!=halo_element_ids.end();++it){
        SURFID[pos++]=*it;
      }
    }
#endif
    return;
  }
  
  void sam_export_halo_fc(int* colgat, int* atosen, int* scater, int* atorec, const int* ncolga, const int* nscate, const int* nprocs, int* pnodes, int* nnodes){
#ifdef HAVE_MPI
    assert(mesh != NULL);
    
    mesh->export_halo(colgat, atosen, scater, atorec, ncolga, nscate, nprocs);

    *pnodes = mesh->node_list.psize();
    *nnodes = mesh->node_list.size();
#endif
    return;
  }

  void sam_export_phalo_fc(int* pcolgat, int* patosen, int* pscater, int* patorec, const int* pncolga, const int* pnscate, const int* nprocs, int* ppnodes, int* pnnodes){
#ifdef HAVE_MPI
    assert(mesh != NULL);
    
    mesh->export_phalo(pcolgat, patosen, pscater, patorec, pncolga, pnscate, nprocs);
    *ppnodes = mesh->MFnode_list.psize();
    *pnnodes = mesh->MFnode_list.size();
#endif
    return;
  }
  
  void sam_cleanup_fc(void){
    if(mesh != NULL){
      delete mesh;
      mesh = NULL;
      
      decomp_opts.clear();
    }
    noddom.clear();
    
    return;
  }

  void sam_add_field_fc(samfloat_t* field_data, int *nnodes){
#ifdef HAVE_MPI
    assert(mesh != NULL);

    assert(mesh->node_list.size() == (size_t)*nnodes);
    for(size_t i=0;i<mesh->node_list.size();i++){
      mesh->node_list[i].append_field(field_data + i, 1);
    }
#endif
    return;
  }
  
  void sam_pop_field_fc(samfloat_t* field_data, int *nnodes){
#ifdef HAVE_MPI
    assert(mesh != NULL);
    
    assert(mesh->node_list.size() == (size_t)*nnodes);
    for (size_t i=0;i<mesh->node_list.size();i++){
      field_data[i] = mesh->node_list[i].pop_field();
    }
#endif
    return;
  }  
  
  void sam_export_node_ownership_fc(int* node_ownership, int* nnodes){
    assert(noddom.size() == *nnodes);
      
    for(int i = 0;i < *nnodes;i++){
      node_ownership[i] = noddom[i];
    }
    
    return;
  }
}
