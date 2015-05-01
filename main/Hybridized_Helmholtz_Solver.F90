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
    type(scalar_field), pointer :: D_rhs, dgified_D, D, &
         dgified_D_rhs, D_exact, D_exact_D_mesh, l_d_rhs
    type(scalar_field), target :: f, D_rhs_projected
    type(vector_field), pointer :: u,X,U_rhs
    type(vector_field) :: U_Local,U_rhs_allocated
    integer :: ierr,dump_no,stat,ele
    character(len = OPTION_PATH_LEN) :: simulation_name
    character(len=PYTHON_FUNC_LEN) :: coriolis
    real :: L2error,Linfty_error,L2projectederror, energy
    logical :: have_d_rhs, have_u_rhs, projection_mode=.false.
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
    D_rhs=>extract_scalar_field(state(1), "LayerThicknessRHS",stat)
    have_d_rhs = (stat==0)
    U_rhs=>extract_vector_field(state(1), "VelocityRHS",stat)
    have_u_rhs = (stat==0)
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

    if(have_d_rhs) then
       call allocate(D_rhs_projected,D%mesh,'D_rhs')
       call remap_field(D_rhs,D_rhs_projected)
       l_d_rhs => d_rhs_projected
    end if

    if(.not.(have_u_rhs).or.(have_d_rhs)) then
       call compute_energy_hybridized(state(1),energy)
       ewrite(1,*) 'ENERGY BEFORE = ', energy
    end if

    dump_no = 0
    call write_state(dump_no,state)

    if(have_d_rhs.or.have_u_rhs) then
       if(.not.have_d_rhs) then
          call solve_hybridized_helmholtz(&
               &state(1),U_rhs=U_rhs,&
               &compute_cartesian=.true.,&
               &check_continuity=.true.,output_dense=.false.,&
               projection=projection_mode)
       else if(have_u_rhs) then
          call solve_hybridized_helmholtz(&
               &state(1),d_rhs=l_d_rhs,U_rhs=U_rhs,&
               &compute_cartesian=.true.,&
               &check_continuity=.true.,output_dense=.false.,&
               projection=projection_mode)          
       else
          call solve_hybridized_helmholtz(&
               &state(1),d_rhs=l_d_rhs,&
               &compute_cartesian=.true.,&
               &check_continuity=.true.,output_dense=.false.,&
               projection=projection_mode)
       end if
    else
       !Timestepping mode.
       ewrite(1,*) 'TIMESTEPPING MODE'
       call solve_hybridized_helmholtz(&
            &state(1),&
            &compute_cartesian=.true.,&
            &check_continuity=.true.,output_dense=.false.,&
            projection=projection_mode)
    end if

    dgified_D => extract_scalar_field(state(1), "LayerThicknessV",stat)
    if(stat==0) then
       call remap_field(D,dgified_D)
    end if
    dgified_D_rhs => extract_scalar_field(state(1), &
         &"LayerThicknessRHSV",stat)
    if(stat==0) then
       call remap_field(D_rhs,dgified_D_rhs)
    end if

    if(projection_mode) then
       call project_local_to_cartesian(X,U_local,U)
       l2error = 0.
       do ele = 1, element_count(U)
          call compute_errors_projection(U,U_rhs,X,ele,l2error)
       end do
       ewrite(1,*) 'l2error for projection', l2error

    else if(have_d_rhs) then
       D_exact => extract_scalar_field(state(1), &
            &"LayerThicknessExact",stat)
       if(stat==0) then
          D_exact_D_mesh => extract_scalar_field(state(1),&
               &"LayerThicknessExactProjected")
          if(stat==0) then
             call remap_field(D_exact,D_exact_D_mesh)
          end if
       end if
       if(stat==0) then
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
       end if
    else
       call compute_energy_hybridized(state(1),energy)
       ewrite(1,*) 'ENERGY AFTER = ', energy
    end if
    call write_state(dump_no,state)

#ifdef HAVE_MPI
    call mpi_finalize(ierr)
    assert(ierr == MPI_SUCCESS)
#endif
  contains
    subroutine compute_errors(D,D_exact,D_exact_D_mesh,X,ele,&
         L2error,Linfty_error,L2projectederror)
      implicit none
      type(scalar_field), intent(in) :: D,D_exact,D_exact_D_mesh
      type(vector_field), intent(in) :: X
      integer, intent(in) :: ele
      real, intent(inout) :: L2error,Linfty_error,L2projectederror
      !
      real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
      real, dimension(ele_ngi(D,ele)) :: D_gi,D_exact_gi,D_exact_D_mesh_gi,&
           &detwei

      call compute_jacobian(X, ele, J=J, detwei=detwei)
      d_gi = ele_val_at_quad(D,ele)
      d_exact_gi = ele_val_at_quad(D_exact,ele)
      d_exact_d_mesh_gi = ele_val_at_quad(D_exact_D_mesh,ele)

      linfty_error = max(linfty_error,maxval(abs(d_gi-d_exact_gi)))
      l2error = l2error + sum((d_gi-d_exact_gi)**2*detwei)
      l2projectederror = l2projectederror + &
           &sum((d_gi-d_exact_d_mesh_gi)**2*detwei)
    end subroutine compute_errors

    subroutine compute_errors_projection(U,U_rhs,X,ele,l2error)
      implicit none
      type(vector_field), intent(in) :: X,U,U_rhs
      integer, intent(in) :: ele
      real, intent(inout) :: L2error
      !
      real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
      real, dimension(ele_ngi(D,ele)) :: detwei
      real, dimension(U%dim, ele_ngi(X,ele)) :: U_quad, U_rhs_quad
      integer :: dim1

      call compute_jacobian(X, ele, J=J, detwei=detwei)
      u_quad = ele_val_at_quad(U,ele)
      u_rhs_quad = ele_val_at_quad(U_rhs,ele)

      do dim1 = 1, U%dim
         l2error = l2error + sum((u_quad(dim1,:)-u_rhs_quad(dim1,:))**2*detwei)
      end do
    end subroutine compute_errors_projection

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

         if(argument=="-p") then
            projection_mode = .true.
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

      write (0,*) "usage: hybridized_helmholtz_solver [-v n] [-p] <options_file>"
      write (0,*) ""
      write (0,*) "-v n sets the verbosity of debugging"
      write (0,*) "-p switches on projection mode (project velocity to div-conforming"
    end subroutine usage

  end program Hybridized_Helmholtz_Solver
