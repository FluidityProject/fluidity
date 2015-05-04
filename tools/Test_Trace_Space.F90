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
  program test_trace_space
    ! Program to test trace space by making some projections between
    ! fields and trace spaces.
    use spud
    use fields
    use state_module
    use FLDebug
    use populate_state_module
    use write_state_module
    use solvers
    use sparse_tools
    use sparsity_patterns_meshes
    use global_parameters, only: option_path_len, python_func_len, current_time, dt
    use memory_diagnostics
    use iso_c_binding
    use mangle_options_tree
#ifdef HAVE_PETSC_MODULES
  use petsc
#endif
  implicit none
#include "petsc_legacy.h"

    ! Interface blocks for the initialisation routines we need to call
    interface
      subroutine set_global_debug_level(n)
        integer, intent(in) :: n
      end subroutine set_global_debug_level

      subroutine python_init
      end subroutine python_init

      subroutine petscinitialize(s, i)
        character(len=*), intent(in) :: s
        integer, intent(out) :: i
      end subroutine petscinitialize
    end interface

    type(state_type), dimension(:), pointer :: state

    integer :: ierr
    character(len = OPTION_PATH_LEN) :: simulation_name

#ifdef HAVE_MPI
    call mpi_init(ierr)
    assert(ierr == MPI_SUCCESS)
#endif

#ifdef HAVE_PETSC
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
#endif
    
    call python_init
    call read_command_line
    call mangle_options_tree_forward
    call set_global_debug_level(1)
    call populate_state(state)
    call get_option('/simulation_name',simulation_name)

    call deallocate_transform_cache

#ifdef HAVE_MEMORY_STATS
    call print_current_memory_stats(0)
#endif

    call test_local_coords(state(1))
    call test_face_local_nodes(state(1))
    call test_trace_values(state(1))
    call test_trace_projection(state(1))

    ewrite(1,*) 'TEST PASSED'

#ifdef HAVE_MPI
    call mpi_finalize(ierr)
    assert(ierr == MPI_SUCCESS)
#endif

  contains

    subroutine test_face_local_nodes(state)
      implicit none
      type(state_type), intent(inout) :: state
      !
      type(scalar_field), pointer :: L
      integer :: ele

      L=>extract_scalar_field(state, "LagrangeMultiplier")

      do ele = 1, ele_count(L)
         call test_face_local_nodes_ele(L,ele)
      end do
      
    end subroutine test_face_local_nodes

    subroutine test_face_local_nodes_ele(L,ele)
      implicit none
      type(scalar_field), intent(in) :: L
      integer, intent(in) :: ele
      !
      integer, dimension(:), pointer :: neigh
      integer :: ni,ele2,face,face2
      
      neigh => ele_neigh(L,ele)
      do ni = 1, size(neigh)
         ele2 = neigh(ni)
         face = ele_face(L,ele,ele2)
         if(ele2>0) then
            face2 = ele_face(L,ele2,ele)
         else
            face2 = -1
         end if
         call test_face_local_nodes_face(L,ele,face,face2)
      end do
    end subroutine test_face_local_nodes_ele

    subroutine test_face_local_nodes_face(L,ele,face,face2)
      implicit none
      type(scalar_field), intent(in) :: L
      integer, intent(in) :: ele,face,face2
      !
      integer, dimension(face_loc(L,face)) :: nods1, nods2, nods3, nods_loc
      integer, dimension(:), pointer :: L_ele

      L_ele => ele_nodes(L,ele)
      nods1 = face_global_nodes(L,face)
      if(face2>0) then
         nods2 = face_global_nodes(L,face2)
      end if
      nods_loc = face_local_nodes(L,face)
      nods3 = L_ele(nods_loc)

      if(face2>0) then
         if(any(nods1/=nods2)) then
            FLAbort('Global node numbers don''t agree on face')
         end if
      end if
      if(any(nods1/=nods3)) then
         FLAbort('Face global nodes doesn''t match global numbers')
      end if
    end subroutine test_face_local_nodes_face

    subroutine test_trace_projection(state)
      implicit none
      !this subroutine checks the values of the local coordinates 
      !for the trace function
      type(state_type), intent(inout) :: state
      !
      type(scalar_field), pointer :: D,L
      type(vector_field), pointer :: X
      type(scalar_field) :: L_projected, L_projected_rhs
      type(csr_sparsity) :: L_mass_sparsity
      type(csr_matrix) :: L_mass_mat
      integer :: i, ele
      D=>extract_scalar_field(state, "LayerThickness")
      L=>extract_scalar_field(state, "LagrangeMultiplier")
      X=>extract_vector_field(state, "Coordinate")
      call allocate(L_projected,L%mesh, "ProjectedLagrangeMultiplier")
      L_projected%option_path = L%option_path
      call allocate(L_projected_rhs,L%mesh,&
           & "ProjectedLagrangeMultiplierRHS")
      call zero(L_projected)
      call zero(L_projected_rhs)

      L_mass_sparsity=get_csr_sparsity_firstorder(state, L%mesh, L%mesh)
      call allocate(L_mass_mat,L_mass_sparsity)
      call zero(L_mass_mat)

      do ele = 1, element_count(L)
         call assemble_trace_projection_ele(ele,D,L_projected_rhs,X,L_mass_mat)
      end do

      call petsc_solve(L_projected, L_mass_mat,&
           &L_projected_rhs)

      if(maxval(abs(L_projected%val-L%val))>1.0e-5) then
         FLExit('Projection to trace space looks funky.')
      end if

      call deallocate(L_mass_mat)
      call deallocate(L_projected) 
     call deallocate(L_projected_rhs)

    end subroutine test_trace_projection

    subroutine assemble_trace_projection_ele(ele,D,L_projected_rhs&
         &,X,L_mass_mat)
      implicit none
      integer, intent(in) :: ele
      type(scalar_field), intent(inout) :: D,L_projected_rhs
      type(vector_field), intent(inout) :: X
      type(csr_matrix), intent(inout) :: L_mass_mat
      !
      integer, dimension(:), pointer :: neigh
      integer :: ni,ele_2,face
      
      neigh => ele_neigh(D,ele)
      do ni = 1, size(neigh)
         ele_2 = neigh(ni)
         if(ele_2<ele) then
            face = ele_face(D,ele,ele_2)
            call assemble_trace_projection_face(face,D,L_projected_rhs,X,&
                 &L_mass_mat)
         end if
      end do

    end subroutine assemble_trace_projection_ele

    subroutine assemble_trace_projection_face(face,D,L_projected_rhs,X&
         &,L_mass_mat)
      implicit none
      integer, intent(in) :: face
      type(scalar_field), intent(inout) :: D,L_projected_rhs
      type(vector_field), intent(inout) :: X
      type(csr_matrix), intent(inout) :: L_mass_mat
      !
      real, dimension(face_loc(D,face),face_loc(D,face)) :: L_mass_mat_face
      real, dimension(face_loc(D,face)) :: L_rhs
      real, dimension(face_ngi(D,face)) :: D_face_quad, detwei
      type(element_type), pointer :: shape
      integer, dimension(face_loc(l_projected_rhs,face)) :: L_face

      l_face = face_global_nodes(L_projected_rhs,face)
      shape => face_shape(L_projected_rhs,face)
      D_face_quad = face_val_at_quad(D,face)
      call transform_facet_to_physical(X, face, &
         &                          detwei_f=detwei)
      L_rhs = shape_rhs(shape,D_face_quad*detwei)
      L_mass_mat_face = shape_shape(shape,shape,detwei)
      
      call addto(L_projected_rhs, l_face, l_rhs)
      call addto(L_mass_mat, l_face, l_face, L_mass_mat_face)

    end subroutine assemble_trace_projection_face

    subroutine test_local_coords(state)
      implicit none
      !this subroutine checks the values of the local coordinates 
      !for the trace function
      type(state_type), intent(inout) :: state
      !
      type(scalar_field), pointer :: D,L
      integer :: i
      D=>extract_scalar_field(state, "LayerThickness")
      L=>extract_scalar_field(state, "LagrangeMultiplier")
      ewrite(2,*) 'L number2count'
      do i = 1, size(L%mesh%shape%numbering%number2count,1)
         ewrite(2,*) L%mesh%shape%numbering%number2count(i,:)
      end do
      ewrite(2,*) 'D local coordinates'
      do i = 1, D%mesh%shape%loc
         ewrite(2,*) i,local_coords(i, D%mesh%shape)
      end do
      ewrite(2,*) 'L local coordinates'
      do i = 1, L%mesh%shape%loc
         ewrite(2,*) i,local_coords(i, L%mesh%shape)
      end do
    end subroutine test_local_coords

    subroutine test_trace_values(state)
      implicit none
      !This subroutine assumes that the layer thickness 
      !is initialised from the same field as the lagrange
      !multiplier and compares values
      type(state_type), intent(inout) :: state
      ! 
      integer :: ele
      type(scalar_field), pointer :: D,L
      D=>extract_scalar_field(state, "LayerThickness")
      L=>extract_scalar_field(state, "LagrangeMultiplier")
      
      do ele = 1, element_count(D)
         print *, 'Testing element ',ele
         call test_trace_values_ele(D,L,ele)
      end do

    end subroutine test_trace_values

    subroutine test_trace_values_ele(D,L,ele)
      implicit none
      type(scalar_field), intent(inout) :: D,L
      integer, intent(in) :: ele
      !
      integer, pointer, dimension(:) :: neigh
      integer :: ni,ele_2,face
      real, pointer, dimension(:) :: D_face, L_face

      print *, 'D_ele', ele_val(D,ele)
      print *, 'L_ele', ele_val(L,ele)
      
      neigh => ele_neigh(D,ele)
      do ni = 1, size(neigh)
         ele_2 = neigh(ni)
         face = ele_face(D,ele,ele_2)
         
         call test_trace_values_face(D,L,face)
      end do
    end subroutine test_trace_values_ele

    subroutine test_trace_values_face(D,L,face)
      implicit none
      type(scalar_field), intent(inout) :: D,L
      integer, intent(in) :: face
      !
      real, dimension(face_loc(D,face)) :: D_face
      real, dimension(face_loc(L,face)) :: L_face

      print *, 'face', face
      
      D_face = face_val(D,face)
      L_face = face_val(L,face)
      
      print *, "D_face", D_face
      print *, "L_face", L_face
      if(any(abs(D_face-L_face)>1.0e-10)) then
         FLExit('Test Trace Values Failed')
      end if
    end subroutine test_trace_values_face

    subroutine read_command_line()
      implicit none
      ! Read the input filename.

      character(len=1024) :: argument
      integer :: status, argn, level

      call set_global_debug_level(0)

      argn=1
      do

         call get_command_argument(argn, value=argument, status=status)
         argn=argn+1

         if (status/=0) then
            call usage
            stop
         end if

         if (argument=="-v") then
            call get_command_argument(argn, value=argument, status=status)
            argn=argn+1

            if (status/=0) then
               call usage
               stop
            end if

            read(argument, "(i1)", err=666) level
            call set_global_debug_level(level)

            ! Go back to pick up the command line.
            cycle
         end if

         exit
      end do

      call load_options(argument)
      if(.not. have_option("/simulation_name")) goto 666

      return

666   call usage
      stop

    end subroutine read_command_line

    subroutine usage
      implicit none

      write (0,*) "usage: test_trace_space [-v n] <options_file>"
      write (0,*) ""
      write (0,*) "-v n sets the verbosity of debugging"
    end subroutine usage

  end program test_trace_space
