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
  program Hybridized_Helmholtz_Solver
    ! Program to solve the Helmholtz problem arising from the implicit
    ! solution of the shallow-water equations.
    ! Uses swml file as input for this reason.
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
    use hybridized_helmholtz
    implicit none
#ifdef HAVE_PETSC
#include "finclude/petsc.h"
#endif

    ! Interface blocks for the initialisation routines we need to call
    interface
      subroutine set_global_debug_level(n)
        integer, intent(in) :: n
      end subroutine set_global_debug_level

      subroutine mpi_init(ierr)
        integer, intent(out) :: ierr
      end subroutine mpi_init

      subroutine mpi_finalize(ierr)
        integer, intent(out) :: ierr
      end subroutine mpi_finalize

      subroutine python_init
      end subroutine python_init

      subroutine petscinitialize(s, i)
        character(len=*), intent(in) :: s
        integer, intent(out) :: i
      end subroutine petscinitialize
    end interface

    type(state_type), dimension(:), pointer :: state
    type(scalar_field), pointer :: D_rhs, dgified_D, D, &
         dgified_D_rhs, D_exact, D_exact_D_mesh
    type(scalar_field) :: f, D_rhs_projected
    type(vector_field), pointer :: u,X,U_rhs
    type(vector_field) :: U_Local
    integer :: ierr,dump_no,stat,ele
    character(len = OPTION_PATH_LEN) :: simulation_name
    character(len=PYTHON_FUNC_LEN) :: coriolis
    real :: L2error,Linfty_error,L2projectederror
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
    D_rhs=>extract_scalar_field(state(1), "LayerThicknessRHS")
    !U_rhs=>extract_vector_field(state(1), "VelocityRHS")
    U=>extract_vector_field(state(1), "Velocity")
    D => extract_scalar_field(state(1), "LayerThickness")
    call allocate(U_local, mesh_dim(U), U%mesh, "LocalVelocity")
    call zero(U_local)
    call insert(state, U_local, "LocalVelocity")
    call deallocate(U_local)

    X=>extract_vector_field(state(1), "Coordinate")
     call allocate(f, X%mesh, "Coriolis")
     call get_option("/physical_parameters/coriolis", coriolis, stat)
    if(stat==0) then
       call set_from_python_function(f, coriolis, X, time=0.0)
    else
       call zero(f)
    end if
    call insert(state, f, "Coriolis")
    call deallocate(f)

    call allocate(D_rhs_projected,D%mesh,'D_rhs')
    call remap_field(D_rhs,D_rhs_projected)
    dump_no = 0
    call write_state(dump_no,state)

    call solve_hybridized_helmholtz(&
         &state(1),D_rhs=D_rhs_projected,&
         &compute_cartesian=.true.,check_continuity=.true.,output_dense=.false.)

    dgified_D => extract_scalar_field(state(1), "LayerThicknessV",stat)
    if(stat==0) then
       call remap_field(D,dgified_D)
    end if
    dgified_D_rhs => extract_scalar_field(state(1), &
         &"LayerThicknessRHSV",stat)
    if(stat==0) then
       call remap_field(D_rhs,dgified_D_rhs)
    end if

    D_exact => extract_scalar_field(state(1), &
         &"LayerThicknessExact")
    D_exact_D_mesh => extract_scalar_field(state(1),&
         &"LayerThicknessExactProjected")
    call remap_field(D_exact,D_exact_D_mesh)

    L2error=0.
    Linfty_error=0.
    L2projectederror=0.
    do ele = 1, ele_count(D)
       call compute_errors(D,D_exact,D_exact_D_mesh,X,ele,&
            L2error,Linfty_error,L2projectederror)
    end do
    l2error = sqrt(l2error)
    l2projectederror = sqrt(l2projectederror)
    ewrite(1,*) 'l_infty error', Linfty_error
    ewrite(1,*) 'l2error', l2error
    ewrite(1,*) 'l2projectederror', l2projectederror
    call write_state(dump_no,state)

#ifdef HAVE_MPI
    call mpi_finalize(ierr)
    assert(ierr == MPI_SUCCESS)
#endif
contains
  subroutine compute_errors(D,D_exact,D_exact_D_mesh,X,ele,&
       L2error,Linfty_error,L2projectederror)
    type(scalar_field), intent(in) :: D,D_exact,D_exact_D_mesh
    type(vector_field), intent(in) :: X
    integer, intent(in) :: ele
    real, intent(inout) :: L2error,Linfty_error,L2projectederror
    !
    real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
    real, dimension(ele_ngi(D,ele)) :: D_gi,D_exact_gi,D_exact_D_mesh_gi,&
         &detwei

    call compute_jacobian(ele_val(X,ele),ele_shape(X,ele), J=J, &
         detwei=detwei)
    d_gi = ele_val_at_quad(D,ele)
    d_exact_gi = ele_val_at_quad(D_exact,ele)
    d_exact_d_mesh_gi = ele_val_at_quad(D_exact_D_mesh,ele)

    linfty_error = max(linfty_error,maxval(abs(d_gi-d_exact_gi)))
    l2error = l2error + sum((d_gi-d_exact_gi)**2*detwei)
    l2projectederror = l2projectederror + &
         &sum((d_gi-d_exact_d_mesh_gi)**2*detwei)
  end subroutine compute_errors

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

      write (0,*) "usage: hybridized_helmholtz_solver [-v n] <options_file>"
      write (0,*) ""
      write (0,*) "-v n sets the verbosity of debugging"
    end subroutine usage

  end program Hybridized_Helmholtz_Solver
