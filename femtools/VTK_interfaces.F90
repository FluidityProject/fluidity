!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineeringp
!    Imperial College London
!
!    C.Pain@Imperial.ac.uk
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

module vtk_interfaces
  ! This module exists purely to provide explicit interfaces to
  ! libvtkfortran. It provides generic interfaces which ensure that the
  ! calls work in both single and double precision.

  use elements
  use fields
  use state_module
  use sparse_tools
  use global_numbering
  use fetools, only: X_,Y_,Z_
  use global_parameters, only : FIELD_NAME_LEN, real_4, real_8
  use mpi_interfaces
  use parallel_tools
  use parallel_fields
  use spud
  
  implicit none

  private

  public :: vtk_write_state, vtk_write_fields, vtk_read_state, &
    vtk_write_surface_mesh
  
  ! Element types from VTK
  integer, public, parameter :: VTK_VERTEX=1
  integer, public, parameter :: VTK_POLY_VERTEX=2
  integer, public, parameter :: VTK_LINE=3
  integer, public, parameter :: VTK_POLY_LINE=4
  integer, public, parameter :: VTK_TRIANGLE=5
  integer, public, parameter :: VTK_TRIANGLE_STRIP=6
  integer, public, parameter :: VTK_POLYGON=7
  integer, public, parameter :: VTK_PIXEL=8
  integer, public, parameter :: VTK_QUAD=9
  integer, public, parameter :: VTK_TETRA=10
  integer, public, parameter :: VTK_VOXEL=11
  integer, public, parameter :: VTK_HEXAHEDRON=12
  integer, public, parameter :: VTK_WEDGE=13
  integer, public, parameter :: VTK_PYRAMID=14
  
  integer, public, parameter :: VTK_QUADRATIC_EDGE=21
  integer, public, parameter :: VTK_QUADRATIC_TRIANGLE=22
  integer, public, parameter :: VTK_QUADRATIC_QUAD=23
  integer, public, parameter :: VTK_QUADRATIC_TETRA=24
  integer, public, parameter :: VTK_QUADRATIC_HEXAHEDRON=25

  interface vtkopen
     SUBROUTINE VTKOPEN(outName, len1, vtkTitle, len2)
       implicit none
       CHARACTER*(*) outName
       INTEGER len1
       CHARACTER vtkTitle
       INTEGER len2
     end SUBROUTINE VTKOPEN
     module procedure vtkopen_f90
  end interface
  
  interface vtkclose
     ! Close the current vtk file.
     subroutine vtkclose()
     end subroutine vtkclose
  end interface

  interface vtkpclose
     ! Close the current vtk file - creates a parallel file.
     subroutine vtkpclose(rank, npartitions)
       implicit none
       integer rank, npartitions
     end subroutine vtkpclose
  end interface

  interface vtkwritemesh
     ! Write mesh information to the current vtk file.
     SUBROUTINE VTKWRITEMESH(NNodes, NElems, x, y, z, enlist, &
          elementTypes, elementSizes)
       use global_parameters, only : real_4
       implicit none
       INTEGER NNodes
       INTEGER NElems
       REAL(real_4) x(*)
       REAL(real_4) y(*)
       REAL(real_4) z(*)
       INTEGER enlist(*)
       INTEGER elementTypes(*)
       INTEGER elementSizes(*)
     end SUBROUTINE VTKWRITEMESH
     SUBROUTINE VTKWRITEMESHD(NNodes, NElems, x, y, z, enlist, &
          elementTypes, elementSizes) 
       use global_parameters, only : real_8
       implicit none
       INTEGER NNodes
       INTEGER NElems
       REAL(real_8) x(*)
       REAL(real_8) y(*)
       REAL(real_8) z(*)
       INTEGER enlist(*)
       INTEGER elementTypes(*)
       INTEGER elementSizes(*)
     end SUBROUTINE VTKWRITEMESHD
  end interface

  interface vtkwritesn
     ! Write a scalar field to the current vtk file.
     SUBROUTINE VTKWRITEISN(vect, name, len)
       implicit none
       INTEGER vect(*)
       CHARACTER*(*) name
       INTEGER len
     end SUBROUTINE VTKWRITEISN
     SUBROUTINE VTKWRITEFSN(vect, name, len)
       use global_parameters, only : real_4
       implicit none
       REAL(real_4) vect(*)
       CHARACTER*(*) name
       INTEGER len
     end SUBROUTINE VTKWRITEFSN
     SUBROUTINE VTKWRITEDSN(vect, name, len)
       use global_parameters, only : real_8
       implicit none
       REAL(real_8) vect(*)
       CHARACTER*(*) name
       INTEGER len
     end SUBROUTINE VTKWRITEDSN
     module procedure vtkwriteisn_f90, vtkwritefsn_f90, vtkwritedsn_f90
  end interface

  interface vtkwritesc
     ! Write a scalar field (cell-based) to the current vtk file.
     SUBROUTINE VTKWRITEISC(vect, name, len)
       implicit none
       INTEGER vect(*)
       CHARACTER*(*) name
       INTEGER len
     end SUBROUTINE VTKWRITEISC
     SUBROUTINE VTKWRITEFSC(vect, name, len)
       use global_parameters, only : real_4
       implicit none
       REAL(real_4) vect(*)
       CHARACTER*(*) name
       INTEGER len
     end SUBROUTINE VTKWRITEFSC
     SUBROUTINE VTKWRITEDSC(vect, name, len)
       use global_parameters, only : real_8
       implicit none
       REAL(real_8) vect(*)
       CHARACTER*(*) name
       INTEGER len
     end SUBROUTINE VTKWRITEDSC
     module procedure vtkwriteisc_f90, vtkwritefsc_f90, vtkwritedsc_f90
  end interface

  interface vtkwritevn
     ! Write a vector field to the current vtk file.
     SUBROUTINE VTKWRITEFVN(vx, vy, vz, name, len)
       use global_parameters, only : real_4
       implicit none
       REAL(real_4) vx(*)
       REAL(real_4) vy(*)
       REAL(real_4) vz(*)
       CHARACTER*(*) name
       INTEGER len
     end SUBROUTINE VTKWRITEFVN
     SUBROUTINE VTKWRITEDVN(vx, vy, vz, name, len)
       use global_parameters, only : real_8
       implicit none
       REAL(real_8) vx(*)
       REAL(real_8) vy(*)
       REAL(real_8) vz(*)
       CHARACTER*(*) name
       INTEGER len
     end SUBROUTINE VTKWRITEDVN
     module procedure vtkwritefvn_f90, vtkwritedvn_f90
  end interface

  interface vtkwritevc
     ! Write a vector field (cell-based) to the current vtk file.
     SUBROUTINE VTKWRITEFVC(vx, vy, vz, name, len)
       use global_parameters, only : real_4
       implicit none
       REAL(real_4) vx(*)
       REAL(real_4) vy(*)
       REAL(real_4) vz(*)
       CHARACTER*(*) name
       INTEGER len
     end SUBROUTINE VTKWRITEFVC
     SUBROUTINE VTKWRITEDVC(vx, vy, vz, name, len)
       use global_parameters, only : real_8
       implicit none
       REAL(real_8) vx(*)
       REAL(real_8) vy(*)
       REAL(real_8) vz(*)
       CHARACTER*(*) name
       INTEGER len
     end SUBROUTINE VTKWRITEDVC
     module procedure vtkwritefvc_f90, vtkwritedvc_f90
  end interface

  interface vtkwritetn
     ! Write a tensor field to the current vtk file.
     SUBROUTINE VTKWRITEFTN(v1, v2, v3, v4, v5, v6, v7, v8, v9, name, len)
       use global_parameters, only : real_4
       implicit none
       REAL(real_4) v1(*)
       REAL(real_4) v2(*)
       REAL(real_4) v3(*)
       REAL(real_4) v4(*)
       REAL(real_4) v5(*)
       REAL(real_4) v6(*)
       REAL(real_4) v7(*)
       REAL(real_4) v8(*)
       REAL(real_4) v9(*)
       CHARACTER*(*) name
       INTEGER len
     end SUBROUTINE VTKWRITEFTN
     SUBROUTINE VTKWRITEDTN(v1, v2, v3, v4, v5, v6, v7, v8, v9, name, len)
       use global_parameters, only : real_8
       implicit none
       REAL(real_8) v1(*)
       REAL(real_8) v2(*)
       REAL(real_8) v3(*)
       REAL(real_8) v4(*)
       REAL(real_8) v5(*)
       REAL(real_8) v6(*)
       REAL(real_8) v7(*)
       REAL(real_8) v8(*)
       REAL(real_8) v9(*)
       CHARACTER*(*) name
       INTEGER len
     end SUBROUTINE VTKWRITEDTN
     module procedure vtkwriteftn_f90, vtkwritedtn_f90
  end interface

  interface vtkwritetc
     ! Write a tensor field (cell-based) to the current vtk file.
     SUBROUTINE VTKWRITEFTC(v1, v2, v3, v4, v5, v6, v7, v8, v9, name, len)
       use global_parameters, only : real_4
       implicit none
       REAL(real_4) v1(*)
       REAL(real_4) v2(*)
       REAL(real_4) v3(*)
       REAL(real_4) v4(*)
       REAL(real_4) v5(*)
       REAL(real_4) v6(*)
       REAL(real_4) v7(*)
       REAL(real_4) v8(*)
       REAL(real_4) v9(*)
       CHARACTER*(*) name
       INTEGER len
     end SUBROUTINE VTKWRITEFTC
     SUBROUTINE VTKWRITEDTC(v1, v2, v3, v4, v5, v6, v7, v8, v9, name, len)
       use global_parameters, only : real_8
       implicit none
       REAL(real_8) v1(*)
       REAL(real_8) v2(*)
       REAL(real_8) v3(*)
       REAL(real_8) v4(*)
       REAL(real_8) v5(*)
       REAL(real_8) v6(*)
       REAL(real_8) v7(*)
       REAL(real_8) v8(*)
       REAL(real_8) v9(*)
       CHARACTER*(*) name
       INTEGER len
     end SUBROUTINE VTKWRITEDTC
     module procedure vtkwriteftc_f90, vtkwritedtc_f90
  end interface

  interface vtkWriteGhostLevels
     subroutine vtkWriteGhostLevels(ghost_levels)
       implicit none
       integer ghost_levels(*)
     end subroutine vtkWriteGhostLevels
  end interface
 
  interface vtksetactivescalars
    subroutine vtksetactivescalars(name, len)
      implicit none
      character*(*) name
      integer len
    end subroutine vtksetactivescalars
    module procedure vtksetactivescalars_f90
  end interface

  interface vtksetactivevectors
    subroutine vtksetactivevectors(name, len)
      implicit none
      character*(*) name
      integer len
    end subroutine vtksetactivevectors
    module procedure vtksetactivevectors_f90
  end interface

  interface vtksetactivetensors
    subroutine vtksetactivetensors(name, len)
      implicit none
      character*(*) name
      integer len
    end subroutine vtksetactivetensors
    module procedure vtksetactivetensors_f90
  end interface

  interface 
       subroutine vtk_read_file(&
          filename, namelen, nnod, nelm, szenls,&
          nfield_components, nprop_components,&
          nfields, nproperties, &
          ndimensions, maxnamelen, &
          x, y, z, &
          field_components, prop_components, &
          fields, properties,&
          enlbas, enlist, field_names, prop_names)
       implicit none
       character*(*) filename
       integer namelen, nnod, nelm, szenls
       integer nfield_components, nprop_components, nfields, nproperties
       integer ndimensions, maxnamelen
       real x(nnod), y(nnod), z(nnod)
       integer field_components(nfields), prop_components(nproperties)
       real fields(nnod,nfields), &
            properties(nelm,nproperties) 
       integer enlbas(nelm+1), enlist(szenls)
       character(len=maxnamelen) field_names(nfields)
       character(len=maxnamelen) prop_names(nproperties)
     end subroutine vtk_read_file
  end interface

  interface
       subroutine vtk_get_sizes( filename, namelen, nnod, nelm, szenls, &
          nfield_components, nprop_components, &
          nfields, nproperties, ndimensions, maxnamelen )
       implicit none
       character*(*) filename
       integer namelen
       integer nnod, nelm, szenls
       integer nfield_components, nprop_components
       integer nfields, nproperties, ndimensions, maxnamelen
     end subroutine vtk_get_sizes
  end interface

contains

  subroutine vtkopen_f90(outName, vtkTitle)
    ! Wrapper routine with nicer interface.
    character(len=*), intent(in) :: outName, vtkTitle

    call vtkopen(outName, len(outName), vtkTitle, len(vtkTitle))

  end subroutine vtkopen_f90
  
  subroutine vtkwriteisn_f90(vect, name)
    ! Wrapper routine with nicer interface.
    integer, intent(in) :: vect(*)
    character(len=*), intent(in) :: name

    call vtkwritesn(vect, name, len(name))

  end subroutine vtkwriteisn_f90

  subroutine vtkwritefsn_f90(vect, name)
    ! Wrapper routine with nicer interface.
    real(real_4), intent(in) :: vect(*)
    character(len=*), intent(in) :: name

    call vtkwritesn(vect, name, len(name))

  end subroutine vtkwritefsn_f90

  subroutine vtkwritedsn_f90(vect, name)
    ! Wrapper routine with nicer interface.
    real(real_8), intent(in) :: vect(*)
    character(len=*), intent(in) :: name

    call vtkwritesn(vect, name, len(name))
  end subroutine vtkwritedsn_f90

  subroutine vtkwriteisc_f90(vect, name)
    ! Wrapper routine with nicer interface.
    integer, intent(in) :: vect(*)
    character(len=*), intent(in) :: name

    call vtkwritesc(vect, name, len(name))

  end subroutine vtkwriteisc_f90

  subroutine vtkwritefsc_f90(vect, name)
    ! Wrapper routine with nicer interface.
    real(real_4), intent(in) :: vect(*)
    character(len=*), intent(in) :: name

    call vtkwritesc(vect, name, len(name))

  end subroutine vtkwritefsc_f90

  subroutine vtkwritedsc_f90(vect, name)
    ! Wrapper routine with nicer interface.
    real(real_8), intent(in) :: vect(*)
    character(len=*), intent(in) :: name

    call vtkwritesc(vect, name, len(name))
  end subroutine vtkwritedsc_f90

  subroutine vtkwritefvn_f90(vx, vy, vz, name)
    REAL(real_4), intent(in) :: vx(*), vy(*), vz(*)
    CHARACTER(len=*) name

    call vtkwritevn(vx, vy, vz, name, len(name))

  end subroutine vtkwritefvn_f90

  subroutine vtkwritedvn_f90(vx, vy, vz, name)
    REAL(real_8), intent(in) :: vx(*), vy(*), vz(*)
    CHARACTER(len=*) name

    call vtkwritevn(vx, vy, vz, name, len(name))

  end subroutine vtkwritedvn_f90  

  subroutine vtkwritefvc_f90(vx, vy, vz, name)
    REAL(real_4), intent(in) :: vx(*), vy(*), vz(*)
    CHARACTER(len=*) name

    call vtkwritevc(vx, vy, vz, name, len(name))

  end subroutine vtkwritefvc_f90

  subroutine vtkwritedvc_f90(vx, vy, vz, name)
    REAL(real_8), intent(in) :: vx(*), vy(*), vz(*)
    CHARACTER(len=*) name

    call vtkwritevc(vx, vy, vz, name, len(name))

  end subroutine vtkwritedvc_f90

  subroutine vtkwriteftn_f90(v1, v2, v3, v4, v5, v6, v7, v8, v9, name)
    REAL(real_4), intent(in) :: v1(*), v2(*), v3(*), v4(*), v5(*), v6(*), v7(*), v8(*), v9(*)
    CHARACTER(len=*) name

    call vtkwritetn(v1, v2, v3, v4, v5, v6, v7, v8, v9, name, len(name))

  end subroutine vtkwriteftn_f90

  subroutine vtkwritedtn_f90(v1, v2, v3, v4, v5, v6, v7, v8, v9, name)
    REAL(real_8), intent(in) :: v1(*), v2(*), v3(*), v4(*), v5(*), v6(*), v7(*), v8(*), v9(*)
    CHARACTER(len=*) name

    call vtkwritetn(v1, v2, v3, v4, v5, v6, v7, v8, v9, name, len(name))

  end subroutine vtkwritedtn_f90  
  
  subroutine vtkwriteftc_f90(v1, v2, v3, v4, v5, v6, v7, v8, v9, name)
    REAL(real_4), intent(in) :: v1(*), v2(*), v3(*), v4(*), v5(*), v6(*), v7(*), v8(*), v9(*)
    CHARACTER(len=*) name

    call vtkwritetc(v1, v2, v3, v4, v5, v6, v7, v8, v9, name, len(name))

  end subroutine vtkwriteftc_f90

  subroutine vtkwritedtc_f90(v1, v2, v3, v4, v5, v6, v7, v8, v9, name)
    REAL(real_8), intent(in) :: v1(*), v2(*), v3(*), v4(*), v5(*), v6(*), v7(*), v8(*), v9(*)
    CHARACTER(len=*) name

    call vtkwritetc(v1, v2, v3, v4, v5, v6, v7, v8, v9, name, len(name))

  end subroutine vtkwritedtc_f90  
  
  subroutine vtksetactivescalars_f90(name)
    character(len=*) name
    
    call vtksetactivescalars(name, len_trim(name))
    
  end subroutine vtksetactivescalars_f90

  subroutine vtksetactivevectors_f90(name)
    character(len=*) name
    
    call vtksetactivevectors(name, len_trim(name))
    
  end subroutine vtksetactivevectors_f90

  subroutine vtksetactivetensors_f90(name)
    character(len=*) name
    
    call vtksetactivetensors(name, len_trim(name))
    
  end subroutine vtksetactivetensors_f90
  
  subroutine vtk_write_state(filename, index, model, state, write_region_ids, stat)
    !!< Write the state variables out to a vtu file. Two different elements
    !!< are supported along with fields corresponding to each of them.
    !!<
    !!< All the fields will be promoted/reduced to the degree of model and
    !!< all elements will be discontinuous (which is required for the
    !!< promotion/reduction to be general).
     implicit none
    character(len=*), intent(in) :: filename !! Base filename with no
    !!< trailing _number.vtu
    integer, intent(in), optional :: index !! Index number of dump for filename.
    character(len=*), intent(in), optional :: model
    type(state_type), dimension(:), intent(in) :: state
    logical, intent(in), optional :: write_region_ids
    integer, intent(out), optional :: stat
    type(mesh_type), pointer :: model_mesh

    ! It is necessary to make local copies of the fields lists because of
    ! the pointer storage mechanism in state
    type(scalar_field), dimension(:), allocatable :: lsfields
    type(vector_field), dimension(:), allocatable :: lvfields
    type(tensor_field), dimension(:), allocatable :: ltfields
    integer :: i, f, counter, size_lsfields, size_lvfields, size_ltfields
    character(len=FIELD_NAME_LEN) :: mesh_name

    if (present(model)) then
      model_mesh => extract_mesh(state(1), model)
    else if (have_option("/io/output_mesh")) then
      ! use the one specified by the options tree:
      call get_option("/io/output_mesh[0]/name", mesh_name)
      model_mesh => extract_mesh(state(1), trim(mesh_name))
    else if (mesh_count(state(1))==1) then
      ! if there's only one mesh, use that:
      model_mesh => extract_mesh(state(1), 1)
    else
      ewrite(-1,*) "In vtk_write_state:"
      FLAbort("Don't know which mesh to use as model.")      
    end if

    size_lsfields = 0
    do i = 1, size(state)
      if (associated(state(i)%scalar_fields)) then
        size_lsfields = size_lsfields + size(state(i)%scalar_fields)
      end if
    end do
    
    allocate(lsfields(size_lsfields))
    counter = 0
    do i = 1, size(state)
      if (associated(state(i)%scalar_fields)) then
        do f = 1, size(state(i)%scalar_fields)
          counter = counter + 1
          lsfields(counter)=state(i)%scalar_fields(f)%ptr       
          if (size(state) > 1) then
            lsfields(counter)%name = trim(state(i)%name)//'::'//trim(lsfields(counter)%name)
          end if
        end do
      end if
    end do

    size_lvfields = 0
    do i = 1, size(state)
      if (associated(state(i)%vector_fields)) then
        size_lvfields = size_lvfields + size(state(i)%vector_fields)
      end if
    end do
    
    allocate(lvfields(size_lvfields))
    counter = 0
    do i = 1, size(state)
      if (associated(state(i)%vector_fields)) then
        do f = 1, size(state(i)%vector_fields)
          counter = counter + 1
          lvfields(counter) = state(i)%vector_fields(f)%ptr
          if (size(state) > 1) then
            lvfields(counter)%name = trim(state(i)%name)//'::'//trim(lvfields(counter)%name)
          end if
        end do
      end if
    end do

    size_ltfields = 0
    do i = 1, size(state)
      if (associated(state(i)%tensor_fields)) then
        size_ltfields = size_ltfields + size(state(i)%tensor_fields)
      end if
    end do
    
    allocate(ltfields(size_ltfields))
    counter = 0
    do i = 1, size(state)
      if (associated(state(i)%tensor_fields)) then
        do f = 1, size(state(i)%tensor_fields)
          counter = counter + 1
          ltfields(counter) = state(i)%tensor_fields(f)%ptr
          if (size(state) > 1) then
            ltfields(counter)%name = trim(state(i)%name)//'::'//trim(ltfields(counter)%name)
          end if
        end do
      end if
    end do

    call vtk_write_fields(filename, index, &
         extract_vector_field(state(1), "Coordinate"), &
         model_mesh,  &
         sfields=lsfields, &
         vfields=lvfields, &
         tfields=ltfields, &
         write_region_ids=write_region_ids, &
         stat=stat)

  end subroutine vtk_write_state

  subroutine vtk_write_fields(filename, index, position, model, sfields,&
       & vfields, tfields, write_region_ids, write_inactive_parts, stat)
    !!< Write the state variables out to a vtu file. Two different elements
    !!< are supported along with fields corresponding to each of them.
    !!<
    !!< All the fields will be promoted/reduced to the degree of model and
    !!< all elements will be discontinuous (which is required for the
    !!< promotion/reduction to be general).
     implicit none

    character(len=*), intent(in) :: filename ! Base filename with no
    ! trailing _number.vtu
    integer, intent(in), optional :: index ! Index number of dump for filename.
    type(vector_field), intent(in) :: position
    type(mesh_type), intent(in) :: model
    type(scalar_field), dimension(:), intent(in), optional :: sfields
    type(vector_field), dimension(:), intent(in), optional :: vfields
    type(tensor_field), dimension(:), intent(in), optional :: tfields
    logical, intent(in), optional :: write_region_ids
    !! If not provided and true, only the local vtu for processes with at least one element are written
    !! The zero element processes are supposed to be trailing in rank. If provided and true all
    !! vtus are written:
    logical, intent(in), optional :: write_inactive_parts
    integer, intent(out), optional :: stat
    
    integer :: NNod, sz_enlist, i, dim, j, k, nparts
    real, dimension(:,:,:), allocatable, target :: t_field_buffer, tensor_values
    real, dimension(:,:), allocatable, target :: v_field_buffer
    real, dimension(:), allocatable, target :: field_buffer
    integer, dimension(:), allocatable, target :: ndglno, ENList, ELsize, ELType
    character(len=FIELD_NAME_LEN) :: dumpnum
    type(mesh_type) :: model_mesh
    type(scalar_field) :: l_model
    type(vector_field) :: v_model(3)
    type(tensor_field) :: t_model
    logical :: dgify_fields ! should we DG-ify the fields -- make them discontinous?
    integer, allocatable, dimension(:)::ghost_levels
    real, allocatable, dimension(:,:) :: tempval
    integer :: lstat
    
    if (present(stat)) stat = 0
    
    dgify_fields = .false.
    if (present(sfields)) then
      do i=1,size(sfields)
        if ( (sfields(i)%mesh%continuity .lt. 0 .and. sfields(i)%mesh%shape%degree /= 0) ) dgify_fields = .true.
      end do
    end if
    if (present(vfields)) then
      do i=1,size(vfields)
        if ( (vfields(i)%mesh%continuity .lt. 0 .and. vfields(i)%mesh%shape%degree /= 0) ) dgify_fields = .true.
      end do
    end if
    if (present(tfields)) then
      do i=1,size(tfields)
        if ( (tfields(i)%mesh%continuity .lt. 0 .and. tfields(i)%mesh%shape%degree /= 0) ) dgify_fields = .true.
      end do
    end if
    
    if (present_and_true(write_inactive_parts)) then
      nparts = getnprocs()
    else
      nparts = get_active_nparts(ele_count(model))
    end if

    if(model%shape%degree /= 0) then
  
      ! Note that the following fails for variable element types.
      sz_enlist = element_count(model)*ele_loc(model, 1)
  
      allocate(ndglno(sz_enlist),ENList(sz_enlist), ELsize(element_count(model)),ELtype(element_count(model)))
  
      if (dgify_fields.and.(continuity(model)>=0)) then
        
        ! Note that the following fails for variable element types.
        NNod=sz_enlist
        
        allocate(field_buffer(NNod), v_field_buffer(NNod, 5), t_field_buffer(NNod, 3, 3))
        v_field_buffer=0.0
        
        call make_global_numbering_DG(NNod, ndglno, element_count(model),&
          & ele_shape(model,1))
      
        ! Discontinuous version of model.
        model_mesh=wrap_mesh(ndglno, ele_shape(model,1), "DGIfiedModelMesh")
        ! Fix bug in gfortran
        model_mesh%shape=ele_shape(model,1)
        ! this mesh is discontinuous
        model_mesh%continuity = -1
        if (associated(model%region_ids)) then
          allocate(model_mesh%region_ids(ele_count(model_mesh)))
          model_mesh%region_ids = model%region_ids
        end if
      else
        
        NNod=node_count(model)
        allocate(field_buffer(NNod), v_field_buffer(NNod, 5), t_field_buffer(NNod, 3, 3))
        v_field_buffer=0.0

        model_mesh = model
        ! Grab an extra reference to make the deallocate at the end safe.
        call incref(model_mesh)
        ndglno = model%ndglno(1:sz_enlist)
      end if
    else
      ! if the model mesh is p/q0 then use the position mesh to output the mesh

      ! Note that the following fails for variable element types.
      sz_enlist = element_count(position%mesh)*ele_loc(position%mesh, 1)
  
      allocate(ndglno(sz_enlist),ENList(sz_enlist), ELsize(element_count(model)),ELtype(element_count(model)))
  
      if(dgify_fields) then
        
        ! Note that the following fails for variable element types.
        NNod=sz_enlist
        
        allocate(field_buffer(NNod), v_field_buffer(NNod, 5), t_field_buffer(NNod, 3, 3))
        v_field_buffer=0.0
        
        call make_global_numbering_DG(NNod, ndglno, element_count(position%mesh),&
          & ele_shape(position%mesh,1))
        ! Discontinuous version of position mesh.
        model_mesh=wrap_mesh(ndglno, ele_shape(position%mesh,1), "")
        ! Fix bug in gfortran
        model_mesh%shape=ele_shape(position%mesh,1)
        ! this mesh is discontinuous
        model_mesh%continuity = -1
        if (associated(model%region_ids)) then
          allocate(model_mesh%region_ids(ele_count(model_mesh)))
          model_mesh%region_ids = model%region_ids
        end if
      else
        
        NNod=node_count(position%mesh)
        allocate(field_buffer(NNod), v_field_buffer(NNod, 5), t_field_buffer(NNod, 3, 3))
        v_field_buffer=0.0
        
        model_mesh = position%mesh
        ! Grab an extra reference to make the deallocate at the end safe.
        call incref(model_mesh)
        ndglno = position%mesh%ndglno
      end if
    end if

    l_model= wrap_scalar_field(model_mesh, field_buffer, "TempVTKModel")
    
    ! vector fields may be of any dimension, so we need a v_model for each case
    ! note that for dim<3, the next columns of v_field_buffer are unused and
    ! always left 0.0. This is needed in the following as vector fields are always
    ! written as 3d (with the unused components being 0)
    v_model(1)=wrap_vector_field(model_mesh,v_field_buffer(:,3), name="TempVTKModel1")
    v_model(2)=wrap_vector_field(model_mesh,v_field_buffer(:,2), &
            v_field_buffer(:,3), name="TempVTKModel2")
    v_model(3)=wrap_vector_field(model_mesh,v_field_buffer(:,1), &
            v_field_buffer(:,2), v_field_buffer(:,3), name="TempVTKModel3")

    t_model=wrap_tensor_field(model_mesh,t_field_buffer, name="TempVTKModel")

    ! Size and type are currently uniform
    ELsize=ele_loc(model_mesh,1)
    ELtype=vtk_element_type(ele_shape(model_mesh,1))

    ENList=fluidity_mesh2vtk_numbering(ndglno, ele_shape(model_mesh,1))


    !----------------------------------------------------------------------
    ! Open the file
    !----------------------------------------------------------------------
    
    if (present(index)) then
       ! Write index number:
       if(nparts > 1) then
          write(dumpnum,"(a,i0,a)") "_",index,".pvtu"
       else
          write(dumpnum,"(a,i0,a)") "_",index,".vtu"
       end if
    else
       ! If no index is provided then assume the filename was complete.
       if(nparts > 1) then
          if (len_trim(filename)<=4) then
             dumpnum=".pvtu"
          else if ((filename(len_trim(filename)-4:len_trim(filename))==".pvtu")) then
             dumpnum=""
          else if (filename(len_trim(filename)-3:len_trim(filename))==".vtu") then
             FLAbort("Parallel vtk write - extension must be .pvtu")
          else
             dumpnum=".pvtu"
          end if
       else
          if (len_trim(filename)<=3) then
             dumpnum=".vtu"
          else if (filename(len_trim(filename)-3:)==".vtu") then
             dumpnum=""
          else
             dumpnum=".vtu"
          end if
       end if
    end if

    if(getprocno() > nparts) then
       return
    end if

    call vtkopen(trim(filename)//trim(dumpnum),trim(filename))

    !----------------------------------------------------------------------
    ! Output the mesh
    !----------------------------------------------------------------------

    ! Remap the position coordinates.
    call remap_field(from_field=position, to_field=v_model(position%dim), stat=lstat)
    ! if this is being called from something other than the main output routines
    ! then these tests can be disabled by passing in the optional stat argument
    ! to vtk_write_fields
    if(lstat==REMAP_ERR_DISCONTINUOUS_CONTINUOUS) then
      if(present(stat)) then
        stat = lstat
      else
        FLAbort("Just remapped from a discontinuous to a continuous field!")
      end if
    else if(lstat==REMAP_ERR_UNPERIODIC_PERIODIC) then
      if(present(stat)) then
        stat = lstat
      else
        ewrite(-1,*) 'While outputting to vtu the coordinates were remapped from'
        ewrite(-1,*) 'a continuous non-periodic to a continuous periodic mesh.'
        ewrite(-1,*) 'This suggests that the output_mesh requested is periodic,'
        ewrite(-1,*) 'which generally produces strange vtus.'
        ewrite(-1,*) 'Best to switch to a non-periodic output_mesh or if you can'
        ewrite(-1,*) "think of a reason why it's a good idea to allow it then delete"
        ewrite(-1,*) "this check (and maybe some that follow) or convert to a warning."
        FLAbort("Just remapped from an unperiodic to a periodic continuous field!")
      end if
    end if
    ! we've just allowed remapping from a higher order to a lower order continuous field

    ! Write the mesh coordinates.
    k=3-position%dim
    call VTKWRITEMESH(node_count(model_mesh), element_count(model_mesh), &
           v_field_buffer(:,k+X_), v_field_buffer(:,k+Y_), v_field_buffer(:,k+Z_)&
           &, ENLIST, ELtype, ELsize)

    !----------------------------------------------------------------------
    ! Output scalar fields
    !----------------------------------------------------------------------    
    
    if (present(sfields)) then
       do i=1,size(sfields)
          
          if (sfields(i)%mesh%shape%degree /= 0) then

            call remap_field(from_field=sfields(i), to_field=l_model, stat=lstat)
            ! if this is being called from something other than the main output routines
            ! then these tests can be disabled by passing in the optional stat argument
            ! to vtk_write_fields
            if(lstat==REMAP_ERR_DISCONTINUOUS_CONTINUOUS) then
              if(present(stat)) then
                stat = lstat
              else
                FLAbort("Just remapped from a discontinuous to a continuous field!")
              end if
            else if(lstat==REMAP_ERR_UNPERIODIC_PERIODIC) then
              if(present(stat)) then
                stat = lstat
              else
                FLAbort("Just remapped from an unperiodic to a periodic continuous field!")
              end if
            end if
            ! we've just allowed remapping from a higher order to a lower order continuous field
            
            call vtkwritesn(l_model%val, trim(sfields(i)%name))

          else

            if(sfields(i)%field_type==FIELD_TYPE_CONSTANT) then
              allocate(tempval(element_count(l_model),1))

              tempval = sfields(i)%val(1)
              call vtkwritesc(tempval(:,1), trim(sfields(i)%name))

              deallocate(tempval)
            else
              call vtkwritesc(sfields(i)%val, trim(sfields(i)%name))
            end if

          end if
          
       end do
       
      ! set first field to be active:
      if (size(sfields)>0) call vtksetactivescalars( sfields(1)%name )
    end if


     !----------------------------------------------------------------------
     ! Output the region ids
     !----------------------------------------------------------------------

     ! You could possibly check for preserving the mesh regions here.
     if (present_and_true(write_region_ids)) then
       if (associated(model_mesh%region_ids)) then
         call vtkwritesc(model_mesh%region_ids, "RegionIds")
       end if
     end if

    !----------------------------------------------------------------------
    ! Output ghost levels
    !----------------------------------------------------------------------
    if(element_halo_count(model_mesh) > 0) then
       allocate(ghost_levels(element_count(model_mesh)))
       ghost_levels = 0
       do i=1, element_count(model_mesh)
         if(.not. element_owned(model, i)) ghost_levels(i) = 1
       end do
       
       call vtkWriteGhostLevels(ghost_levels)
    end if

    
    !----------------------------------------------------------------------
    ! Output vector fields
    !----------------------------------------------------------------------    
    
    if (present(vfields)) then
       do i=1,size(vfields)
          if(trim(vfields(i)%name)=="Coordinate") then
             cycle
          end if

          if(vfields(i)%mesh%shape%degree /= 0) then

            call remap_field(from_field=vfields(i), to_field=v_model(vfields(i)%dim), stat=lstat)
            ! if this is being called from something other than the main output routines
            ! then these tests can be disabled by passing in the optional stat argument
            ! to vtk_write_fields
            if(lstat==REMAP_ERR_DISCONTINUOUS_CONTINUOUS) then
              if(present(stat)) then
                stat = lstat
              else
                FLAbort("Just remapped from a discontinuous to a continuous field!")
              end if
            else if(lstat==REMAP_ERR_UNPERIODIC_PERIODIC) then
              if(present(stat)) then
                stat = lstat
              else
                FLAbort("Just remapped from an unperiodic to a periodic continuous field!")
              end if
            end if
            ! we've just allowed remapping from a higher order to a lower order continuous field
          
            k=3-vfields(i)%dim
            call vtkwritevn(&
                v_field_buffer(:,k+X_), v_field_buffer(:,k+Y_), &
                v_field_buffer(:,k+Z_), &
                trim(vfields(i)%name))

          else

            allocate(tempval(element_count(model_mesh),3))

            tempval = 0.0
            if(vfields(i)%field_type==FIELD_TYPE_CONSTANT) then
              do j = 1, vfields(i)%dim
                tempval(:,j) = vfields(i)%val(j)%ptr(1)
              end do
            else
              do j = 1, vfields(i)%dim
                tempval(:,j) = vfields(i)%val(j)%ptr(:)
              end do
            end if

            call vtkwritevc(&
                tempval(:,X_), tempval(:,Y_), &
                tempval(:,Z_), trim(vfields(i)%name))

            deallocate(tempval)

          end if

       end do

      ! set first field to be active:
      do i=1,size(vfields)
         if(trim(vfields(i)%name)=="Coordinate") then
             cycle
         end if
         call vtksetactivevectors( vfields(i)%name )
         exit
      end do
    end if


    !----------------------------------------------------------------------
    ! Output tensor fields
    !----------------------------------------------------------------------    
    
    if (present(tfields)) then

       do i=1,size(tfields)
          dim = tfields(i)%dim
          
          if(tfields(i)%mesh%shape%degree /= 0) then

            call remap_field(from_field=tfields(i), to_field=t_model, stat=lstat)
            ! if this is being called from something other than the main output routines
            ! then these tests can be disabled by passing in the optional stat argument
            ! to vtk_write_fields
            if(lstat==REMAP_ERR_DISCONTINUOUS_CONTINUOUS) then
              if(present(stat)) then
                stat = lstat
              else
                FLAbort("Just remapped from a discontinuous to a continuous field!")
              end if
            else if(lstat==REMAP_ERR_UNPERIODIC_PERIODIC) then
              if(present(stat)) then
                stat = lstat
              else
                FLAbort("Just remapped from an unperiodic to a periodic continuous field!")
              end if
            end if
            ! we've just allowed remapping from a higher order to a lower order continuous field
            
            allocate(tensor_values(node_count(t_model), 3, 3))
            tensor_values=0.0
            do j=1,dim
              do k=1,dim
                tensor_values(:, j, k) = t_model%val(j, k, :)
              end do
            end do
  
            call vtkwritetn(tensor_values(:, 1, 1), &
                            tensor_values(:, 1, 2), &
                            tensor_values(:, 1, 3), &
                            tensor_values(:, 2, 1), &
                            tensor_values(:, 2, 2), &
                            tensor_values(:, 2, 3), &
                            tensor_values(:, 3, 1), &
                            tensor_values(:, 3, 2), &
                            tensor_values(:, 3, 3), &
                            trim(tfields(i)%name))
            deallocate(tensor_values)

        else

            allocate(tensor_values(element_count(t_model), 3, 3))
            tensor_values=0.0
            if(tfields(i)%field_type==FIELD_TYPE_CONSTANT) then
              do j=1,dim
                do k=1,dim
                  tensor_values(:, j, k) = tfields(i)%val(j, k, 1)
                end do
              end do
            else
              do j=1,dim
                do k=1,dim
                  tensor_values(:, j, k) = tfields(i)%val(j, k, :)
                end do
              end do
            end if

            call vtkwritetc(tensor_values(:, 1, 1), &
                            tensor_values(:, 1, 2), &
                            tensor_values(:, 1, 3), &
                            tensor_values(:, 2, 1), &
                            tensor_values(:, 2, 2), &
                            tensor_values(:, 2, 3), &
                            tensor_values(:, 3, 1), &
                            tensor_values(:, 3, 2), &
                            tensor_values(:, 3, 3), &
                            trim(tfields(i)%name))

            deallocate(tensor_values)

        end if

       end do
       
      ! set first field to be active:
      if (size(tfields)>0) call vtksetactivetensors( tfields(1)%name )
      
    end if


    !----------------------------------------------------------------------
    ! Close the file
    !----------------------------------------------------------------------
    call deallocate(l_model)
    do i=1, 3
      call deallocate(v_model(i))
    end do
    call deallocate(t_model)
    call deallocate(model_mesh)

    if(nparts > 1) then
       call vtkpclose(getrank(), nparts)
    else
       call vtkclose()
    end if
    
  end subroutine vtk_write_fields

  function fluidity_mesh2vtk_numbering(ndglno, element) result (renumber)
    type(element_type), intent(in) :: element
    integer, dimension(:), intent(in) :: ndglno
    integer, dimension(size(ndglno)) :: renumber

    integer, dimension(element%loc) :: ele_num
    integer :: i, nloc

    ele_num=vtk2fluidity_ordering(element)

    nloc=element%loc

    forall (i=1:size(ndglno)/nloc)
       renumber((i-1)*nloc+1:i*nloc)=ndglno((i-1)*nloc+ele_num)
    end forall

  end function fluidity_mesh2vtk_numbering

  function vtk_mesh2fluidity_numbering(vtk_ndglno, element) result (fl_ndglno)
    type(element_type), intent(in) :: element
    integer, dimension(:), intent(in) :: vtk_ndglno
    integer, dimension(size(vtk_ndglno)) :: fl_ndglno

    integer, dimension(element%loc) :: ele_num
    integer :: i, nloc

    ele_num=vtk2fluidity_ordering(element)

    nloc=element%loc

    forall (i=1:size(vtk_ndglno)/nloc)
       fl_ndglno((i-1)*nloc+ele_num)=vtk_ndglno((i-1)*nloc+1:i*nloc)
    end forall

  end function vtk_mesh2fluidity_numbering
  
  function vtk_element_type(element) result (type)
    ! Return the vtk element type corresponding to element.
    ! return 0 if no match is found.
    integer :: type
    type(element_type), intent(in) :: element

    type=0

    select case (element%dim)
    case (1)
       ! Interval elements.
       select case (element%numbering%degree)
       case (0)
          type=VTK_VERTEX
       case (1)
          type=VTK_LINE
       case(2)
          type=VTK_QUADRATIC_EDGE
       case default
          ewrite(0,*) "Polynomial degree: ", element%numbering%degree
          FLAbort("Unsupported polynomial degree for vtk.")
       end select
    case(2)
       select case(element%numbering%vertices)
       case (3)
          select case (element%numbering%degree)
          case (0)
             type=VTK_VERTEX
          case(1)
             type=VTK_TRIANGLE
          case(2)
             type=VTK_QUADRATIC_TRIANGLE
          case default
             ewrite(0,*) "Polynomial degree: ", element%numbering%degree
             FLAbort("Unsupported polynomial degree for vtk.")
          end select
       case (4)
          select case (element%numbering%degree)
          case (0)
             type=VTK_VERTEX
          case(1)
             type=VTK_QUAD
          case(2)
             type=VTK_QUADRATIC_QUAD
          case default
             ewrite(0,*) "Polynomial degree: ", element%numbering%degree
             FLAbort("Unsupported polynomial degree for vtk.")
          end select
       case default
          ewrite(0,*) "Dimension: ", element%dim
          ewrite(0,*) "Vertices: ", element%numbering%vertices
          FLAbort("Unsupported element type for vtk.")
       end select
    case(3)
       select case(element%numbering%vertices)
       case (4)
          select case (element%numbering%degree)
          case (0)
             type=VTK_VERTEX
          case(1)
             type=VTK_TETRA
          case(2)
             type=VTK_QUADRATIC_TETRA
          case default
             ewrite(0,*) "Polynomial degree: ", element%numbering%degree
             FLAbort("Unsupported polynomial degree for vtk.")
          end select
       case (8)
          select case (element%numbering%degree)
          case (0)
             type=VTK_VERTEX
          case(1)
             type=VTK_HEXAHEDRON
          case(2)
             type=VTK_QUADRATIC_HEXAHEDRON
          case default
             ewrite(0,*) "Polynomial degree: ", element%numbering%degree
             FLAbort("Unsupported polynomial degree for vtk.")
          end select
        case default
          ewrite(0,*) "Dimension: ", element%dim
          ewrite(0,*) "Vertices: ", element%numbering%vertices
          FLAbort("Unsupported element type for vtk.")
       end select
    case default
       ewrite(0,*) "Dimension: ", element%dim
       FLAbort("Unsupported dimension for vtk.")
    end select

  end function vtk_element_type

  function vtk2fluidity_ordering(element) result (order)
    ! Return the One True Element Numbering for element relative to the VTK ordering.
    !
    ! Note that the one true element numbering does not contain information
    ! on chirality so transformed elements may have the oposite chirality
    ! to that expected by VTK.
    type(element_type), intent(in) :: element
    integer, dimension(element%loc) :: order
    
    integer :: type

    type=vtk_element_type(element)

    order=0

    select case(type)
    case(VTK_VERTEX)
       order=(/1/)
    case(VTK_LINE)
       order=(/1,2/)
    case(VTK_QUADRATIC_EDGE)
       order=(/1,3,2/)
    case(VTK_TRIANGLE)
       order=(/1,2,3/)
    case(VTK_QUADRATIC_TRIANGLE)
       order=(/1,3,6,2,5,4/)
    case(VTK_QUAD)
       ! this is already reordered inside libvtkfortran :(
       order=(/1,2,3,4/)
    case(VTK_TETRA)
       order=(/1,2,3,4/)
    case(VTK_QUADRATIC_TETRA)
       order=(/1,3,6,10,2,5,4,7,8,9/)
    case(VTK_HEXAHEDRON)
       ! this is already reordered inside libvtkfortran :(
       order=(/1,2,3,4,5,6,7,8/)
    ! NOTE: quadratic quads and hexes are not supported as
    ! vtk quadratic quads/hexes are only quadratic along the edges
    ! i.e. there are no internal nodes.
    case default
       ewrite(0,*) "VTK element type: ", type
       FLAbort("Unsupported element type")
    end select

  end function vtk2fluidity_ordering

  subroutine vtk_read_state(filename, state, quad_degree)
    !!<  This routine uses the vtkmeshio operations
    !!<  to extract mesh and field information from a VTU file.
    character(len=*), intent(in) :: filename
    type(state_type), intent(inout) :: state
    integer, intent(in), optional :: quad_degree

    type(quadrature_type) :: quad
    type(element_type) :: shape
    type(mesh_type) :: mesh, p0_mesh
    type(vector_field) :: position_field

    integer :: i

    integer :: nodes, elements, dim, sz_enlist
    integer :: nfields, nprops, nfield_components, nprop_components
    integer :: degree, maxnamelen
    real, allocatable :: X(:), Y(:), Z(:)
    real, dimension(:,:), allocatable :: fields, properties
    integer, dimension(:), allocatable :: field_components, prop_components
    integer, allocatable :: ENLBAS(:), NDGLNO(:)
    character(len=FIELD_NAME_LEN), allocatable :: field_names(:), prop_names(:)

    integer :: nloc, quaddegree, nvertices

    call nullify(state)

    ! needed for fgetvtksizes, to tell it to work out the dimension
    dim = 0
    
    call vtk_get_sizes(trim(filename), len_trim(filename), nodes, elements,&
         & sz_enlist, nfield_components, nprop_components, &
         & nfields, nprops, dim, maxnamelen)

    nloc = sz_enlist/elements
    ! set quadrature to max available
    select case(dim)
    case(3)
      select case (nloc)
      case (4, 10) ! tets
        quaddegree = 8
      case (8, 27) ! hexes
        quaddegree = 7
      case default
        FLAbort("Unknown element type!")
      end select
    case(2)
      select case (nloc)
      case (3, 6) ! triangles
        quaddegree = 8
      case (4, 9) ! quads
        quaddegree = 9
      case default
        FLAbort("Unknown element type!")
      end select
    case(0)
      ewrite(-1, *) "For vtu filename: " // trim(filename)
      FLExit("vtu not found")
    case default
      ewrite(-1, *) "Dimension: ", dim
      FLAbort("Only 2 or 3 dimensions implemented currently.")
    end select
    if (present(quad_degree)) then
       quaddegree=quad_degree
    end if
    
    allocate(X(nodes), Y(nodes), Z(nodes))
    allocate(FIELDS(nodes, nfield_components), PROPERTIES(elements, nprop_components))
    allocate(field_components(nfields), prop_components(nprops))
    allocate(ENLBAS(elements+1), NDGLNO(sz_enlist))
    allocate(field_names(nfields), prop_names(nprops))

    do i=1, nfields
      field_names(i) = ' '
    end do
    do i=1, nprops
      prop_names(i) = ' '
    end do

    call vtk_read_file(trim(filename), len_trim(filename), &
         & nodes, elements, sz_enlist, &
         & nfield_components, nprop_components, &
         & nfields, nprops, dim, FIELD_NAME_LEN, &
         & X, Y, Z, &
         & field_components, prop_components, &
         & FIELDS, PROPERTIES, &
         & ENLBAS, NDGLNO, &
         & field_names, prop_names)

    if (nloc == 10 .and. dim==3) then
       nvertices=4 ! quadratic tets
       degree=2
    else if (nloc == 6 .and. dim==2) then
       nvertices=3 ! quadratic triangles
       degree=2
    else if (nloc == 27 .and. dim==3) then
       nvertices=8 ! quadratic hexes
       degree=2
    else if (nloc == 9 .and. dim==2) then
       nvertices=4 ! quadratic quads
       degree=2
    else
       ! linear:
       nvertices=nloc
       degree=1
    end if
    
    quad = make_quadrature(vertices=nvertices, dim=dim, degree=quaddegree)
    shape = make_element_shape(vertices=nvertices, dim=dim, degree=degree, quad=quad)
    call allocate(mesh, nodes, elements, shape, name="Mesh")
    mesh%ndglno=vtk_mesh2fluidity_numbering(ndglno, shape)
    call deallocate(shape)
    ! heuristic check for discontinous meshes
    if (nloc*elements==nodes .and. elements > 1) then
      mesh%continuity=-1
    end if
    call insert(state, mesh, "Mesh")

    call allocate(position_field, dim, mesh, name="Coordinate")
    call set_all(position_field, 1, x)
    if (dim>1) then
      call set_all(position_field, 2, y)
    end if
    if (dim>2) then
      call set_all(position_field, 3, z)
    end if

    call insert(state, position_field, "Coordinate")
    
    if (nprops>0) then
      ! cell-wise data is stored in the arrays properties(:,1:nprops)
      ! this is returned as fields on a p0 mesh
      shape = make_element_shape(vertices=nvertices, dim=dim, degree=0, quad=quad)
      p0_mesh=make_mesh(mesh, shape=shape, continuity=-1, name="P0Mesh")
      call deallocate(shape)
      call insert(state, p0_mesh, name="P0Mesh")
    end if
    
    ! insert point-wise fields
    call vtk_insert_fields_in_state(state, &
      mesh, field_components, fields, field_names, dim)
      
    if (nprops>0) then
      ! insert cell-wise fields
      call vtk_insert_fields_in_state(state, &
        p0_mesh, prop_components, properties, prop_names, dim)
    end if

    deallocate(enlbas, ndglno, field_names, prop_names)
    deallocate(properties, fields)
    deallocate(field_components, prop_components)
    deallocate(x, y, z)
    call deallocate(quad)
    call deallocate(mesh)
    if (nprops>0) call deallocate(p0_mesh)
    call deallocate(position_field)

  end subroutine vtk_read_state
    
  subroutine vtk_insert_fields_in_state(state, &
    mesh, components, fields, names, ndim)
  ! insert the fields returned by vtk_read_file in state
    type(state_type), intent(inout):: state
    integer, dimension(:), intent(in):: components
    type(mesh_type), intent(inout):: mesh
    real, dimension(:,:):: fields
    character(len=*), dimension(:), intent(in):: names
    integer, intent(in):: ndim
    
    type(tensor_field):: tfield
    type(vector_field):: vfield
    type(scalar_field):: sfield
    integer:: i, j, k, component, ndim2
    
    component = 1
    do i=1, size(names)
      
      if (components(i)==9 .or. components(i)==4) then
        
        if (components(i)==9) then
          ndim2=3
        else
          ndim2=2
        end if
        ! Let's make a tensor field, see?
        call allocate(tfield, mesh, names(i))
        call zero(tfield)
        do j=1, ndim2
          do k=1, ndim2
            if (j<=ndim .and. k<=ndim) then
              call set_all(tfield, dim1=j, dim2=k, &
                   val=fields(:, component))
            end if
            component = component+1
          end do
        end do
        call insert(state, tfield, names(i))
        call deallocate(tfield)
        
      else if (components(i)==2 .or. components(i)==3) then
        ! Let's make a vector field.
        call allocate(vfield, ndim, mesh, NAMES(i))
        call zero(vfield)
        do j=1, components(i)
          if (j<=ndim) then
            call set_all(vfield, dim=j, val=fields(:, component))
          end if
          component = component+1
        end do
        call insert(state, vfield, NAMES(i))
        call deallocate(vfield)
        
      else if (components(i)==1) then
        ! a scalar field
        call allocate(sfield, mesh, names(i))
        call set_all(sfield, fields(:, component))
        call insert(state, sfield, names(i))
        component = component+1
        call deallocate(sfield)
        
      else
      
        ewrite(-1,*) "In vtk_read_state ***"
        ewrite(-1,*) "Field ", trim(names(i)), " has ", components(i), " components."
        FLAbort("Don't know what to do with that number of components")
        
      end if
      
    end do
    
  end subroutine vtk_insert_fields_in_state
    
  subroutine vtk_write_surface_mesh(filename, index, position)
    character(len=*), intent(in):: filename
    integer, intent(in), optional:: index
    type(vector_field), intent(in), target:: position
      
    type(vector_field):: surface_position
    type(scalar_field), dimension(:), allocatable:: sfields
    type(mesh_type), pointer:: mesh
    type(mesh_type):: pwc_mesh
    integer, dimension(:), allocatable:: surface_element_list
    integer:: i
    
    mesh => position%mesh    
    
    assert( has_faces(mesh) )
    pwc_mesh = piecewise_constant_mesh(mesh%faces%surface_mesh, "PWCSurfaceMesh")
    if (associated(mesh%faces%coplanar_ids)) then
      allocate( sfields(1:2) )
    else
      allocate( sfields(1) )
    end if
    
    call allocate( sfields(1), pwc_mesh, name="BoundaryIDs")
    call set(sfields(1), (/ ( i, i=1,node_count(pwc_mesh)) /), &
      float(mesh%faces%boundary_ids))
    
    if (associated(mesh%faces%coplanar_ids)) then
      call allocate( sfields(2), pwc_mesh, name="CoplanarIDs")
      call set(sfields(2), (/ ( i, i=1,node_count(pwc_mesh)) /), &
        float(mesh%faces%coplanar_ids))
    end if
    call deallocate(pwc_mesh)
    
    allocate(surface_element_list(1:surface_element_count(position)))
    call allocate(surface_position, position%dim, mesh%faces%surface_mesh, "SurfacePositions")
    do i=1, surface_element_count(position)
      surface_element_list(i)=i
    end do
    call remap_field_to_surface(position, surface_position, surface_element_list)
        
    ! some domains may have no surface elements, so we need write_inactive_parts=.true.
    call vtk_write_fields(filename, index=index, position=surface_position, &
      model=mesh%faces%surface_mesh, sfields=sfields, write_inactive_parts=.true.)
      
    call deallocate(sfields(1))
    if (associated(mesh%faces%coplanar_ids)) then
      call deallocate(sfields(2))
    end if
    call deallocate(surface_position)
    
    deallocate(sfields, surface_element_list)    
    
  end subroutine vtk_write_surface_mesh

end module vtk_interfaces
