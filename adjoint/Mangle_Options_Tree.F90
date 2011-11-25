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
module mangle_options_tree
  use spud
  use global_parameters, only: OPTION_PATH_LEN
  use field_options
  use fldebug
  use futils
  implicit none

  private
  public :: mangle_options_tree_forward, mangle_options_tree_adjoint, adjoint_field_path

  contains

  subroutine mangle_options_tree_forward
    ! Change the options dictionary to make it suitable
    ! for the forward model.
    integer :: field, state, functional, nfields, nfields_to_delete
    character(len=OPTION_PATH_LEN) :: path, field_name
    character(len=OPTION_PATH_LEN), dimension(:), allocatable :: fields_to_delete
    integer :: stat
    integer :: j

    do state=0,option_count("/material_phase")-1
      nfields = option_count("/material_phase[" // int2str(state) // "]/scalar_field")
      allocate(fields_to_delete(nfields))

      ! We need to delete fields that are marked as only living in the adjoint.
      nfields_to_delete = 0

      do field=0,nfields-1
        path = "/material_phase[" // int2str(state) // "]/scalar_field["//int2str(field)//"]"
        call get_option(trim(path)//"/name", field_name)

        if (have_option(trim(complete_field_path(trim(path))) // '/adjoint_storage/exists_in_adjoint')) then
          ! We can't actually delete it now, because if we have scalar_field[0,1,2]
          ! and delete scalar_field[1],
          ! now the numbering will be scalar_field[0,1]
          ! and the loop will go all wrong.
          nfields_to_delete = nfields_to_delete + 1
          fields_to_delete(nfields_to_delete) = field_name
        end if
      end do

      do field=1,nfields_to_delete
        path = "/material_phase[" // int2str(state) // "]/scalar_field::" // trim(fields_to_delete(field))
        call delete_option(trim(path))
      end do

      deallocate(fields_to_delete)

      nfields = option_count("/material_phase[" // int2str(state) // "]/vector_field")
      allocate(fields_to_delete(nfields))

      ! We need to delete fields that are marked as only living in the adjoint.
      nfields_to_delete = 0

      do field=0,nfields-1
        path = "/material_phase[" // int2str(state) // "]/vector_field["//int2str(field)//"]"
        call get_option(trim(path)//"/name", field_name)

        if (have_option(trim(complete_field_path(trim(path))) // '/adjoint_storage/exists_in_adjoint')) then
          ! We can't actually delete it now, because if we have vector_field[0,1,2]
          ! and delete vector_field[1],
          ! now the numbering will be vector_field[0,1]
          ! and the loop will go all wrong.
          nfields_to_delete = nfields_to_delete + 1
          fields_to_delete(nfields_to_delete) = field_name
        end if

        ! We should also check if any field has a /no_initial_condition -- it must not
        do j=0,option_count(trim(complete_field_path(trim(path))) // "/initial_condition/no_initial_condition")-1
          if (have_option(trim(complete_field_path(trim(path))) // "/initial_condition[" // int2str(j) // &
                         & "]/no_initial_condition")) then
            FLExit("Users must not specify /no_initial_condition for themselves (path=" // trim(path))
          endif
        end do
      end do

      do field=1,nfields_to_delete
        path = "/material_phase[" // int2str(state) // "]/vector_field::" // trim(fields_to_delete(field))
        call delete_option(trim(path))
      end do

      deallocate(fields_to_delete)

      nfields = option_count("/material_phase[" // int2str(state) // "]/tensor_field")
      allocate(fields_to_delete(nfields))

      ! We need to delete fields that are marked as only living in the adjoint.
      nfields_to_delete = 0

      do field=0,nfields-1
        path = "/material_phase[" // int2str(state) // "]/tensor_field["//int2str(field)//"]"
        call get_option(trim(path)//"/name", field_name)

        if (have_option(trim(complete_field_path(trim(path))) // '/adjoint_storage/exists_in_adjoint')) then
          ! We can't actually delete it now, because if we have tensor_field[0,1,2]
          ! and delete tensor_field[1],
          ! now the numbering will be tensor_field[0,1]
          ! and the loop will go all wrong.
          nfields_to_delete = nfields_to_delete + 1
          fields_to_delete(nfields_to_delete) = field_name
        end if
      end do

      do field=1,nfields_to_delete
        path = "/material_phase[" // int2str(state) // "]/tensor_field::" // trim(fields_to_delete(field))
        call delete_option(trim(path))
      end do

      deallocate(fields_to_delete)
    end do

  end subroutine mangle_options_tree_forward

  subroutine mangle_options_tree_adjoint
    ! Change the options dictionary to make it suitable
    ! for the adjoint model

    integer :: field, state, nfields, nfields_to_delete, nfields_to_move
    character(len=OPTION_PATH_LEN) :: path, field_name, simulation_name
    character(len=OPTION_PATH_LEN), dimension(:), allocatable :: fields_to_delete, fields_to_move
    integer :: stat
    real :: finish_time, current_time, dt

    do state=0,option_count("/material_phase")-1
      nfields = option_count("/material_phase[" // int2str(state) // "]/scalar_field")
      allocate(fields_to_delete(nfields))

      ! We need to delete fields that are marked as only living in the forward run.
      nfields_to_delete = 0

      do field=0,nfields-1
        path = "/material_phase[" // int2str(state) // "]/scalar_field["//int2str(field)//"]"
        call get_option(trim(path)//"/name", field_name)

        if (have_option(trim(complete_field_path(trim(path))) // '/adjoint_storage/exists_in_forward')) then
          ! We can't actually delete it now, because if we have scalar_field[0,1,2]
          ! and delete scalar_field[1],
          ! now the numbering will be scalar_field[0,1]
          ! and the loop will go all wrong.
          nfields_to_delete = nfields_to_delete + 1
          fields_to_delete(nfields_to_delete) = field_name
        end if
      end do

      do field=1,nfields_to_delete
        path = "/material_phase[" // int2str(state) // "]/scalar_field::" // trim(fields_to_delete(field))
        call delete_option(trim(path))
      end do

      deallocate(fields_to_delete)

      nfields = option_count("/material_phase[" // int2str(state) // "]/vector_field")
      allocate(fields_to_delete(nfields))

      ! We need to delete fields that are marked as only living in the adjoint.
      nfields_to_delete = 0

      do field=0,nfields-1
        path = "/material_phase[" // int2str(state) // "]/vector_field["//int2str(field)//"]"
        call get_option(trim(path)//"/name", field_name)

        if (have_option(trim(complete_field_path(trim(path))) // '/adjoint_storage/exists_in_forward')) then
          ! We can't actually delete it now, because if we have vector_field[0,1,2]
          ! and delete vector_field[1],
          ! now the numbering will be vector_field[0,1]
          ! and the loop will go all wrong.
          nfields_to_delete = nfields_to_delete + 1
          fields_to_delete(nfields_to_delete) = field_name
        end if
      end do

      do field=1,nfields_to_delete
        path = "/material_phase[" // int2str(state) // "]/vector_field::" // trim(fields_to_delete(field))
        call delete_option(trim(path))
      end do

      deallocate(fields_to_delete)

      nfields = option_count("/material_phase[" // int2str(state) // "]/tensor_field")
      allocate(fields_to_delete(nfields))

      ! We need to delete fields that are marked as only living in the adjoint.
      nfields_to_delete = 0

      do field=0,nfields-1
        path = "/material_phase[" // int2str(state) // "]/tensor_field["//int2str(field)//"]"
        call get_option(trim(path)//"/name", field_name)

        if (have_option(trim(complete_field_path(trim(path))) // '/adjoint_storage/exists_in_forward')) then
          ! We can't actually delete it now, because if we have tensor_field[0,1,2]
          ! and delete tensor_field[1],
          ! now the numbering will be tensor_field[0,1]
          ! and the loop will go all wrong.
          nfields_to_delete = nfields_to_delete + 1
          fields_to_delete(nfields_to_delete) = field_name
        end if
      end do

      do field=1,nfields_to_delete
        path = "/material_phase[" // int2str(state) // "]/tensor_field::" // trim(fields_to_delete(field))
        call delete_option(trim(path))
      end do

      deallocate(fields_to_delete)
    end do


    do state=0,option_count("/material_phase")-1
      nfields = option_count("/material_phase[" // int2str(state) // "]/scalar_field")
      allocate(fields_to_move(nfields))
      nfields_to_move = 0

      ! Now we need to go through and change the names of all prognostic fields
      ! to AdjointOldName.
      do field=0,nfields-1
        path = "/material_phase[" // int2str(state) // "]/scalar_field["//int2str(field)//"]"
        call get_option(trim(path) // "/name", field_name)
        if (have_option(trim(path) // "/prognostic")) then
          ! We can't actually move it now, for exactly the same reason as above
          nfields_to_move = nfields_to_move + 1
          fields_to_move(nfields_to_move) = field_name
        end if
      end do

      do field=1,nfields_to_move
        path = "/material_phase[" // int2str(state) // "]/scalar_field::"
        call move_option(trim(path) // trim(fields_to_move(field)), trim(path) // "Adjoint" // trim(fields_to_move(field)), stat=stat)
        assert(stat == SPUD_NO_ERROR)
        ! Delete the sources specified for the forward model -- if desired
        if (have_option(trim(path) // "Adjoint" // trim(fields_to_move(field)) // "/prognostic/scalar_field::Source/")) then
          if (have_option(trim(complete_field_path(trim(path) // "Adjoint" // trim(fields_to_move(field)) // "/prognostic/scalar_field::Source/")) &
                        & // "adjoint_storage/exists_in_forward")) then
            call delete_option(trim(path) // "Adjoint" // trim(fields_to_move(field)) // "/prognostic/scalar_field::Source", stat=stat)
          end if
        end if
        ! And delete any initial conditions -- we will also specify these
        if (have_option(trim(path) // "Adjoint" // trim(fields_to_move(field)) // "/prognostic/initial_condition")) then
          call delete_option(trim(path) // "Adjoint" // trim(fields_to_move(field)) // "/prognostic/initial_condition", stat=stat)
          call add_option(trim(path) // "Adjoint" // trim(fields_to_move(field)) // "/prognostic/initial_condition/no_initial_condition", stat=stat)
        endif
      end do

      deallocate(fields_to_move)

      nfields = option_count("/material_phase[" // int2str(state) // "]/vector_field")
      allocate(fields_to_move(nfields))
      nfields_to_move = 0

      ! Now we need to go through and change the names of all prognostic fields
      ! to AdjointOldName.
      do field=0,nfields-1
        path = "/material_phase[" // int2str(state) // "]/vector_field["//int2str(field)//"]"
        call get_option(trim(path) // "/name", field_name)
        if (have_option(trim(path) // "/prognostic")) then
          ! We can't actually move it now, for exactly the same reason as above
          nfields_to_move = nfields_to_move + 1
          fields_to_move(nfields_to_move) = field_name
        end if
      end do

      do field=1,nfields_to_move
        path = "/material_phase[" // int2str(state) // "]/vector_field::"
        call move_option(trim(path) // trim(fields_to_move(field)), trim(path) // "Adjoint" // trim(fields_to_move(field)), stat=stat)
        assert(stat == SPUD_NO_ERROR)
        ! Delete any sources specified for the forward model -- we deal with these ourselves
        if (have_option(trim(path) // "Adjoint" // trim(fields_to_move(field)) // "/prognostic/vector_field::Source/")) then
          if (have_option(trim(complete_field_path(trim(path) // "Adjoint" // trim(fields_to_move(field)) // "/prognostic/vector_field::Source/")) &
                        & // "adjoint_storage/exists_in_forward")) then
            call delete_option(trim(path) // "Adjoint" // trim(fields_to_move(field)) // "/prognostic/vector_field::Source", stat=stat)
          end if
        end if
        ! And delete any initial conditions -- we will also specify these
        if (have_option(trim(path) // "Adjoint" // trim(fields_to_move(field)) // "/prognostic/initial_condition")) then
          call delete_option(trim(path) // "Adjoint" // trim(fields_to_move(field)) // "/prognostic/initial_condition", stat=stat)
          call add_option(trim(path) // "Adjoint" // trim(fields_to_move(field)) // "/prognostic/initial_condition/no_initial_condition", stat=stat)
        endif
      end do

      deallocate(fields_to_move)

      nfields = option_count("/material_phase[" // int2str(state) // "]/tensor_field")
      allocate(fields_to_move(nfields))
      nfields_to_move = 0

      ! Now we need to go through and change the names of all prognostic fields
      ! to AdjointOldName.
      do field=0,nfields-1
        path = "/material_phase[" // int2str(state) // "]/tensor_field["//int2str(field)//"]"
        call get_option(trim(path) // "/name", field_name)
        if (have_option(trim(path) // "/prognostic")) then
          ! We can't actually move it now, for exactly the same reason as above
          nfields_to_move = nfields_to_move + 1
          fields_to_move(nfields_to_move) = field_name
        end if
      end do

      do field=1,nfields_to_move
        path = "/material_phase[" // int2str(state) // "]/tensor_field::"
        call move_option(trim(path) // trim(fields_to_move(field)), trim(path) // "Adjoint" // trim(fields_to_move(field)), stat=stat)
        assert(stat == SPUD_NO_ERROR)
        ! Delete any sources specified for the forward model -- we deal with these ourselves
        if (have_option(trim(path) // "Adjoint" // trim(fields_to_move(field)) // "/prognostic/tensor_field::Source/")) then
          if (have_option(trim(complete_field_path(trim(path) // "Adjoint" // trim(fields_to_move(field)) // "/prognostic/tensor_field::Source/")) &
                        & // "adjoint_storage/exists_in_forward")) then
            call delete_option(trim(path) // "Adjoint" // trim(fields_to_move(field)) // "/prognostic/tensor_field::Source", stat=stat)
          end if
        end if
        ! And delete any initial conditions -- we will also specify these
        if (have_option(trim(path) // "Adjoint" // trim(fields_to_move(field)) // "/prognostic/initial_condition")) then
          call delete_option(trim(path) // "Adjoint" // trim(fields_to_move(field)) // "/prognostic/initial_condition", stat=stat)
          call add_option(trim(path) // "Adjoint" // trim(fields_to_move(field)) // "/prognostic/initial_condition/no_initial_condition", stat=stat)
        endif
      end do

      deallocate(fields_to_move)
    end do

    ! And the simulation name
    call get_option("/simulation_name", simulation_name)
    call set_option("/simulation_name", trim(simulation_name) // "_adjoint", stat=stat)

    ! And the timestepping information
    call get_option("/timestepping/current_time", current_time)
    call get_option("/timestepping/finish_time", finish_time)
    call get_option("/timestepping/timestep", dt)

    call set_option("/timestepping/current_time", finish_time)
    call set_option("/timestepping/finish_time", current_time)
    call set_option("/timestepping/timestep", -dt)

  end subroutine mangle_options_tree_adjoint

  function adjoint_field_path(fwd_path) result(adj_path)
    character(len=*), intent(in) :: fwd_path
    character(len=OPTION_PATH_LEN) :: adj_path

    integer :: j, k
    integer :: no_fwd_slash
    integer :: no_colons

    adj_path = ""
    no_fwd_slash = 0
    j = 1
    k = 1
    ! First mode: find the second slash in /material_phase::whatever/scalar_field::whatever
    do while(j <= OPTION_PATH_LEN .and. k <= len_trim(fwd_path) .and. no_fwd_slash /= 2)
      adj_path(j:j) = fwd_path(k:k)
      if (fwd_path(k:k) == "/") then
        no_fwd_slash = no_fwd_slash + 1
      end if
      j = j + 1
      k = k + 1
    end do

    no_colons = 0

    ! Now that we are here, we need to move to the next :: in between this slash and the next / (if any)
    do while (j <= OPTION_PATH_LEN .and. k <= len_trim(fwd_path) .and. no_fwd_slash == 2 .and. no_colons /= 2)
      adj_path(j:j) = fwd_path(k:k)
      if (fwd_path(k:k) == "/") then
        no_fwd_slash = no_fwd_slash + 1
      else if (fwd_path(k:k) == ':') then
        no_colons = no_colons + 1
      end if
      j = j + 1
      k = k + 1
    end do

    ! Now that we are here, we have either hit the next /, and can just run on, or
    ! if no_colons == 2, we need to insert an Adjoint
    if (no_colons == 2) then
      adj_path(j:j+6) = "Adjoint"
      j = j + 7
    end if

    do while(k <= len_trim(fwd_path))
      adj_path(j:j) = fwd_path(k:k)
      j = j + 1
      k = k + 1
    end do
  end function adjoint_field_path
end module mangle_options_tree
