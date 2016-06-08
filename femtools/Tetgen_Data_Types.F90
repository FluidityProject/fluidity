!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    amcgsoftware@imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#include "fdebug.h"

module tetgen_data_types

  use iso_c_binding

  implicit none

  private

  type, bind(c) :: mesh_data
     integer(c_int) :: nnodes, nelements, nfacets, nholes, nfaces, nattributes
     integer(c_int) :: lnode_ids, lregion_ids, lface_ids
     type(c_ptr)    :: nodes 
     type(c_ptr)    :: node_ids
     type(c_ptr)    :: ndglno
     type(c_ptr)    :: region_ids
     type(c_ptr)    :: region_attributes
     type(c_ptr)    :: facets
     type(c_ptr)    :: face_ids
     type(c_ptr)    :: holes
     type(c_ptr)    :: element_adjacency
     type(c_ptr)    :: faces
     type(c_ptr)    :: faces_adjacency
  end type mesh_data

  public mesh_data

end module tetgen_data_types
