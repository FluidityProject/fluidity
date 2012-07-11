!    Copyright (C) 2009 Imperial College London and others.
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
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
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
  
  module bubble_tools
    use fields
    use fields_allocates
    use vector_tools
    use sparse_matrices_fields
    use solvers
    use quadrature
    use sparsity_patterns_meshes
    use element_numbering
    use state_module
    use shape_functions
    use quadrature
    use shape_functions
    use fldebug_parameters
    use ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use spud
    use vtk_interfaces
    use adjacency_lists
    use reference_counting
    implicit none

    private
    public get_lumped_mass_p2b,&
         & project_to_p2b_lumped, get_p2b_lumped_mass_quadrature, &
         & bubble_field_to_vtk

#include "../femtools/Reference_count_interface_mesh_type.F90"

  contains

    !! This is special cased because it is only used for mass lumping
    subroutine get_p2b_lumped_mass_quadrature(quad)
      type(quadrature_type), intent(inout) :: quad
      !
      real :: wv = 1.0/20., we = 2.0/15., wg = 9./20.

      call allocate(quad, 3, 7, 3)
      quad%dim = 2
      quad%degree = 3
      quad%weight = 0.5*(/wv,we,wv,we,we,wv,wg/)
      quad%l(:,1) = (/1.0,0.5,0.0,0.5,0.0,0.,1./3./)
      quad%l(:,2) = (/0.0,0.5,1.0,0.0,0.5,0.,1./3./)
      quad%l(:,3) = (/0.0,0.0,0.0,0.5,0.5,1.,1./3./)
      quad%family = 3!FAMILY_SIMPSONS

    end subroutine get_p2b_lumped_mass_quadrature

    !! This is special cased because of the special properties of the 
    !! p2b lumped mass (it is still 3rd order, and positive).
    subroutine get_lumped_mass_p2b(state,mass,p2b_field)
      type(csr_matrix), intent(inout) :: mass
      type(state_type), intent(inout) :: state
      type(scalar_field), intent(in) :: p2b_field
      !
      integer :: ele
      type(vector_field), pointer :: X
      real, dimension(7) :: N_vals
      type(quadrature_type) :: Simpsons_quad
      type(element_type) :: p2b_lumped_shape, X_lumped_shape

      X=>extract_vector_field(state, "Coordinate")

      if(p2b_field%mesh%shape%dim.ne.2) then
         FLAbort('Only works for 2d meshes')
      end if
      if(p2b_field%mesh%shape%loc.ne.7) then
         FLAbort('Expected p2 bubble mesh')
      end if

      call get_p2b_lumped_mass_quadrature(Simpsons_quad)
      p2b_lumped_shape = make_element_shape_from_element(p2b_field%mesh%shape, &
           quad=Simpsons_quad)
      X_lumped_shape = make_element_shape_from_element(X%mesh%shape, &
           quad=Simpsons_quad)

      call zero(mass)

      do ele = 1, ele_count(X)
         call get_lumped_mass_p2b_ele(mass,p2b_field,p2b_lumped_shape,&
              X,X_lumped_shape,&
              ele)
      end do
      
      call deallocate(p2b_lumped_shape)
      call deallocate(X_lumped_shape)
      call deallocate(Simpsons_quad)

    end subroutine get_lumped_mass_p2b

    subroutine get_lumped_mass_p2b_ele(mass,p2b_field,p2b_lumped_shape,&
         X,X_lumped_shape,ele)
      type(csr_matrix), intent(inout) :: mass
      type(scalar_field), intent(in) :: p2b_field
      type(vector_field), intent(in) :: X
      type(element_type), intent(in) :: p2b_lumped_shape,X_lumped_shape
      integer, intent(in) :: ele
      !
      real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
      real, dimension(ele_ngi(X,ele)) :: detwei
      real, dimension(X_lumped_shape%ngi) :: detwei_l
      real, dimension(ele_loc(p2b_field,ele),ele_loc(p2b_field,ele))&
           :: l_mass_mat
      real :: Area
      integer :: loc
      real, dimension(p2b_lumped_shape%ngi) :: weight_vals
      call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J, detwei)
      Area = sum(detwei)
      call compute_jacobian(ele_val(X,ele), X_lumped_shape, J, detwei_l)
      assert(abs(Area-sum(detwei_l))<1.0e-10)
      l_mass_mat = shape_shape(p2b_lumped_shape,p2b_lumped_shape,&
           detwei_l)
      call addto(mass,ele_nodes(p2b_field,ele),&
           ele_nodes(p2b_field,ele),l_mass_mat)

    end subroutine get_lumped_mass_p2b_ele

    !! This is special cased because of the special properties of the 
    !! p2b lumped mass (it is still 3rd order, and positivity preserving)
    subroutine project_to_p2b_lumped(state,field,field_projected)
      type(state_type), intent(inout) :: state
      type(scalar_field), intent(in) :: field
      type(scalar_field), intent(inout) :: field_projected
      !
      type(scalar_field) :: rhs
      type(vector_field), pointer :: X
      type(csr_sparsity), pointer :: mass_sparsity
      type(csr_matrix) :: lumped_mass
      integer :: ele
      type(quadrature_type) :: Simpsons_quad
      type(element_type) :: X_lumped_shape,&
           & field_projected_lumped_shape, field_lumped_shape

      if(field_projected%mesh%shape%dim.ne.2) then
         FLAbort('Only works for 2d meshes')
      end if
      if(field_projected%mesh%shape%loc.ne.7) then
         FLAbort('Expected p2 bubble mesh')
      end if
      X=>extract_vector_field(state, "Coordinate")

      call get_p2b_lumped_mass_quadrature(Simpsons_quad)
      field_projected_lumped_shape = make_element_shape_from_element( &
           field_projected%mesh%shape,quad=Simpsons_quad)
      field_lumped_shape = make_element_shape_from_element( &
           field%mesh%shape,quad=Simpsons_quad)
      X_lumped_shape = make_element_shape_from_element(X%mesh%shape, &
           quad=Simpsons_quad)

      mass_sparsity => get_csr_sparsity_firstorder(state, &
           field_projected%mesh, field_projected%mesh)
      call allocate(lumped_mass,mass_sparsity)
      call zero(lumped_mass)
      call allocate(rhs,field_projected%mesh,"ProjectionRHS")
      call zero(rhs)

      do ele = 1, ele_count(X)
         call project_to_p2b_lumped_ele(rhs,lumped_mass,&
              X,X_lumped_shape,field_lumped_shape,&
              field_projected_lumped_shape,field,ele)
      end do
      call petsc_solve(field_projected,lumped_mass,rhs)
      ewrite (2,*) maxval(abs(field_projected%val))
      call deallocate(rhs)
      call deallocate(lumped_mass)
    end subroutine project_to_p2b_lumped

    subroutine project_to_p2b_lumped_ele(rhs,lumped_mass,&
         X,X_lumped_shape,field_lumped_shape,&
         field_projected_lumped_shape,field,ele)
      type(vector_field), intent(in) :: X
      type(scalar_field), intent(in) :: field
      type(scalar_field), intent(inout) :: rhs
      type(csr_matrix), intent(inout) :: lumped_mass
      integer, intent(in) :: ele
      type(element_type), intent(in) :: X_lumped_shape,&
           field_projected_lumped_shape, field_lumped_shape
      !
      real, dimension(ele_loc(rhs,ele)) :: l_rhs,field_gi      
      real, dimension(ele_loc(rhs,ele),ele_loc(rhs,ele)) :: &
           & l_mass
      real, dimension(ele_loc(field,ele)) :: field_vals
      integer :: loc
      real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
      real, dimension(X_lumped_shape%ngi) :: detwei

      call compute_jacobian(ele_val(X,ele), X_lumped_shape, J, detwei)

      !Values of field at field DOFs
      field_vals = ele_val(field,ele)
      assert(size(field_gi)==size(detwei))
      field_gi = matmul(transpose(field_lumped_shape%n),field_vals)
      l_rhs = shape_rhs(field_projected_lumped_shape,&
           field_gi*detwei)
      l_mass = shape_shape(field_projected_lumped_shape,&
           field_projected_lumped_shape,detwei)
      call addto(rhs,ele_nodes(rhs,ele),l_rhs)
      call addto(lumped_mass,ele_nodes(rhs,ele),&
           ele_nodes(rhs,ele),l_mass)
    end subroutine project_to_p2b_lumped_ele

    subroutine bubble_field_to_vtk(state,sfields,filename,dump_num)
      type(state_type), intent(inout) :: state
      type(scalar_field), dimension(:), intent(in) :: sfields
      character(len=*), intent(in) :: filename ! Base filename
      integer, intent(in) :: dump_num
      !
      type(mesh_type), pointer :: lQuadMesh
      type(mesh_type), target :: QuadMesh
      type(vector_field) :: QuadX
      type(scalar_field), dimension(size(sfields)) :: QuadSF
      type(vector_field), pointer :: X
      integer :: stat, i

      ewrite(1,*) 'subroutine bubble_field_to_vtk(state,sfield,filename,dump_num)'
      do i = 1, size(sfields)
         if(sfields(i)%mesh%shape%numbering%type.ne.ELEMENT_BUBBLE) then
            FLExit('Only works for bubble functions.')
         end if
      end do
      X => extract_vector_field(state,"Coordinate")

      !Construct a quadrilateral mesh from the triangular one.
      lQuadMesh => extract_mesh(state,"QuadMesh",stat)
      if(stat/=0) then
         call quadmesh_from_trimesh(QuadMesh,X%mesh,state)
         lQuadMesh => QuadMesh
      end if
      !Construct a coordinate field on the quadrilateral mesh
      call allocate(QuadX,X%dim,lQuadMesh,"QuadCoordinate")
      call zero(QuadX)
      call quadX_from_bubbleX(QuadX,X)
      !Remap the bubble field to Q1 on the quad mesh
      do i = 1, size(sfields)
         call allocate(QuadSF(i),lQuadMesh,"Quad"//trim(sfields(i)%name))
         call quadSF_from_bubbleSF(QuadSF(i),sfields(i))
      end do
      call vtk_write_fields(filename, dump_num, quadX, lQuadMesh, &
           QuadSF)
      call deallocate(QuadX)
      do i = 1, size(sfields)
         call deallocate(QuadSF(i))
      end do
      
      ewrite(1,*) 'END subroutine bubble_field_to_vtk(state,sfield,filename,dump_num)'
      contains 

        subroutine quadmesh_from_trimesh(QuadMesh,Xmesh,state)
          type(mesh_type), intent(inout) :: QuadMesh
          type(mesh_type), intent(in) :: Xmesh
          type(state_type), intent(inout) :: state
          !
          type(csr_sparsity) :: NNList, EEList
          type(csr_matrix) :: edge_matrix
          type(quadrature_type) :: quad
          integer :: tri_edges, tri_nodes, tri_elements
          integer :: ele, ni, nodecount, node, qele, qele_loc,ele2
          integer, pointer, dimension(:) :: neigh, X_ele

          tri_elements = Xmesh%elements
          quadmesh%elements = 3*tri_elements
          allocate(quadmesh%ndglno(quadmesh%elements*4))
          quadmesh%ndglno=-666
          call MakeLists_Mesh(Xmesh,NNList=NNlist,EEList=EEList)
          call allocate(edge_matrix,EEList,type=CSR_INTEGER)
          call zero(edge_matrix)
          tri_edges=size(NNList%colm)/2
          tri_nodes=Xmesh%nodes
          quadmesh%nodes = 3*tri_elements+tri_edges+Xmesh%nodes

          nodecount = Xmesh%nodes
          do ele = 1, ele_count(Xmesh)             
             !! Element numbering in the three quads
             !! ele_q = 3*(ele-1)+vnode where vnode is the local
             !! node number of the node at the vertex
             !!Local numbering in the three quads
             !!triangle vertex node first, then edge node
             !!on edge joining to vertex with next number in modulo arithmetic
             !!then triangle centre, then other edge

             !First do the vertex nodes (same numbering as Xmesh)
             X_ele => ele_nodes(Xmesh,ele)
             assert(size(X_ele)==3)
             do qele_loc = 1,3
                qele = 3*(ele-1)+qele_loc
                quadmesh%ndglno(4*(qele-1)+1) = X_ele(qele_loc)
             end do
             
             !Then the edge nodes (only do if ele<ele2 or ele2<0)
             neigh => ele_neigh(Xmesh,ele)
             assert(size(neigh)==3)
             do ni = 1, size(neigh)
                ele2 = neigh(ni)
                if(ele2.le.0 .or. ele<ele2) then
                   nodecount = nodecount + 1
                   node = nodecount
                   if(ele2>0) then
                      call set(edge_matrix,ele,ele2,node)
                      call set(edge_matrix,ele2,ele,node)
                   end if
                else
                   node = val(edge_matrix,ele,ele2)
                end if
                assert(node>0)
                select case(ni)
                case (1) 
                   qele = 3*(ele-1)+2
                   quadmesh%ndglno(4*(qele-1)+2)=node
                   qele = 3*(ele-1)+3
                   quadmesh%ndglno(4*(qele-1)+3)=node
                case (2)
                   qele = 3*(ele-1)+3
                   quadmesh%ndglno(4*(qele-1)+2)=node
                   qele = 3*(ele-1)+1
                   quadmesh%ndglno(4*(qele-1)+3)=node
                case (3)
                   qele = 3*(ele-1)+1
                   quadmesh%ndglno(4*(qele-1)+2)=node
                   qele = 3*(ele-1)+2
                   quadmesh%ndglno(4*(qele-1)+3)=node
                end select
             end do

             !Then the triangle centre nodes
             do qele = 3*(ele-1)+1,3*(ele-1)+3
                nodecount = nodecount + 1
                quadmesh%ndglno(4*(qele-1)+4)=nodecount
             end do
          end do
          assert(all(quadmesh%ndglno>0))
          assert(nodecount==quadmesh%nodes)

          quadmesh%wrapped=.false.
          quad = make_quadrature(4,2,degree=Xmesh%shape%quadrature%degree,&
               family=Xmesh%shape%quadrature%family)
          quadmesh%shape = make_element_shape(4,2,1,quad)
          quadmesh%name = "QuadMesh"
          quadmesh%continuity=0
          quadmesh%periodic=.false.

          allocate(quadmesh%adj_lists)
          nullify(quadmesh%region_ids)
          nullify(quadmesh%subdomain_mesh)
          nullify(quadmesh%refcount) ! Hack for gfortran component initialisation
          !                         bug.

          call addref(quadmesh)
          
          !need to insert into state and deallocate
          call insert(state,quadmesh,"QuadMesh")
          call deallocate(quadmesh)
          call deallocate(EEList)
          call deallocate(NNList)
        end subroutine quadmesh_from_trimesh

        subroutine quadX_from_bubbleX(QuadX,X)
          type(vector_field), intent(inout) :: QuadX
          type(vector_field), intent(in) :: X
          !
          real, dimension(X%dim,ele_loc(X,1)) :: X_vals
          real, dimension(X%dim,ele_loc(QuadX,1)) :: QX_vals
          integer :: qele,dim1,ele

          do ele = 1, ele_count(X)
             X_vals = ele_val(X,ele)

             !! Quadrilateral #1
             do dim1 = 1, X%dim
                QX_vals(dim1,:) = &
                     &(/X_vals(dim1,1),0.5*(X_vals(dim1,1)+X_vals(dim1,2)),&
                     &0.5*(X_vals(dim1,1)+X_vals(dim1,3)),&
                     &sum(X_vals(dim1,:))/3/)
             end do
             qele = (ele-1)*3+1
             call set(QuadX,ele_nodes(QuadX,qele),QX_vals)
             !! Quadrilateral #2
             do dim1 = 1, X%dim
                QX_vals(dim1,:) = &
                     &(/X_vals(dim1,2),0.5*(X_vals(dim1,2)+X_vals(dim1,3)),&
                     &0.5*(X_vals(dim1,1)+X_vals(dim1,2)),&
                     &sum(X_vals(dim1,:))/3/)
             end do
             qele = (ele-1)*3+2
             call set(QuadX,ele_nodes(QuadX,qele),QX_vals)                
             !! Quadrilateral #3
             do dim1 = 1, X%dim
                QX_vals(dim1,:) = &
                     &(/X_vals(dim1,3),0.5*(X_vals(dim1,3)+X_vals(dim1,1)),&
                     &0.5*(X_vals(dim1,2)+X_vals(dim1,3)),&
                     &sum(X_vals(dim1,:))/3/)
             end do
             qele = (ele-1)*3+3
             call set(QuadX,ele_nodes(QuadX,qele),QX_vals)                
          end do

        end subroutine quadX_from_bubbleX

        subroutine quadSF_from_bubbleSF(QuadSF,sfield)
          type(scalar_field), intent(inout) :: QuadSF
          type(scalar_field), intent(in) :: sfield
          
          !
          real, dimension(ele_loc(sfield,1)) :: s_vals
          real, dimension(ele_loc(QuadSF,1)) :: qs_vals
          integer :: qele,ele
          real, dimension(7) :: N_vals
          
          !Basis functions evaluated at bubble node.
          N_vals = eval_shape(ele_shape(sfield,1), (/1.0/3.0,1.0/3.0,1.0/3.0/))

          do ele = 1, ele_count(sfield)
             s_vals = ele_val(sfield,ele)
             assert(size(s_vals)==7)
             s_vals(7) = s_vals(7)/N_vals(7)
             s_vals(7) = s_vals(7) + sum(s_vals(1:6)*N_vals(1:6))
             !! Quadrilateral #1
             QS_vals(:) = &
                  &(/s_vals(1),s_vals(2),s_vals(4),s_vals(7)/)
             qele = (ele-1)*3+1
             call set(QuadSF,ele_nodes(QuadSF,qele),QS_vals)
             !! Quadrilateral #2
             QS_vals(:) = &
                  &(/s_vals(3),s_vals(5),s_vals(2),s_vals(7)/)
             qele = (ele-1)*3+2
             call set(QuadSF,ele_nodes(QuadSF,qele),QS_vals)
             !! Quadrilateral #3
             QS_vals(:) = &
                  &(/s_vals(6),s_vals(4),s_vals(5),s_vals(7)/)
             qele = (ele-1)*3+3
             call set(QuadSF,ele_nodes(QuadSF,qele),QS_vals)
          end do
        end subroutine quadSF_from_bubbleSF
      end subroutine bubble_field_to_vtk

#include "../femtools/Reference_count_mesh_type.F90"
      
  end module bubble_tools
