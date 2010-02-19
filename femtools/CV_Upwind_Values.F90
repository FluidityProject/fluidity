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
module cv_upwind_values
  !!< Module containing general tools for discretising Control Volume problems.
  use quadrature
  use elements
  use spud
  use fields
  use sparse_tools
  use state_module
  use fldebug
  use cv_shape_functions, only: make_cvbdy_element_shape
  use cv_faces
  use cvtools, only: complete_cv_field_path
  use cv_options
  use vector_tools, only: norm2, cross_product, scalar_triple_product
  use transform_elements, only: transform_cvsurf_facet_to_physical
  use field_derivatives, only: grad
  use global_parameters, only: OPTION_PATH_LEN
  use boundary_conditions, only: get_periodic_boundary_condition, &
                                 get_entire_boundary_condition, &
                                 get_boundary_condition_nodes

  implicit none

  ! critical distance to project out of element
  real, private, parameter :: c_distance=0.001
  real, private, parameter :: tolerance=tiny(0.0)

  interface find_upwind_values
    module procedure find_upwind_values_single_state, find_upwind_values_multiple_states
  end interface

  private
  public :: need_upwind_values, &
            find_upwind_values, &
            calculate_boundary_normals, &
            couple_upwind_values

contains

  logical function need_upwind_values(option_path)
    ! This function checks whether a field with option_path requires the calculation of
    ! upwind values its control volume spatial discretisation.

    character(len=*), intent(in) :: option_path


    need_upwind_values = ((have_option(trim(complete_cv_field_path(option_path))//&
                          "/face_value[0]/limit_face_value")).or.&
                          (have_option(trim(complete_cv_field_path(option_path))//&
                          "/face_value::HyperC")).or.&
                          (have_option(trim(complete_cv_field_path(option_path))//&
                          "/face_value::UltraC")).or.&
                          (have_option(trim(complete_cv_field_path(option_path))//&
                          "/face_value::PotentialUltraC")))

  end function need_upwind_values
  
  subroutine find_upwind_values_single_state(state, x_field, field, upwind_values, &
                                old_field, old_upwind_values, &
                                defer_deletion, option_path)
    ! This subroutine wraps the various methods for calculating upwind values.
    ! It returns a csr matrix (which must be allocated before) with the upwind
    ! values for each node pair.

    ! bucket full of fields
    type(state_type), intent(inout) :: state
    ! the coordinates on a mesh similar (i.e. not the same when periodic) to the field
    type(vector_field), intent(inout) :: x_field
    ! the field and its previous timelevel that we're interested in
    type(scalar_field), intent(inout) :: field, old_field
    ! the matrices of upwind and old upwind values
    type(csr_matrix), intent(inout) :: upwind_values, old_upwind_values
    ! do we want to temporarily insert the upwind elements and quadrature matrices
    ! into state so that other calls from this subroutine can use them?
    logical, optional :: defer_deletion
    ! for back compatibility pass in an option_path in case field is
    ! locally wrapped or allocated
    character(len=*), optional :: option_path
  
    type(state_type), dimension(1) :: states
  
    states = (/state/)
    call find_upwind_values(states, x_field, field, upwind_values, &
                                old_field, old_upwind_values, &
                                defer_deletion=defer_deletion, option_path=option_path)
    state = states(1)
    
  end subroutine find_upwind_values_single_state
  
  subroutine find_upwind_values_multiple_states(state, x_field, field, upwind_values, &
                                old_field, old_upwind_values, &
                                defer_deletion, option_path)
    ! This subroutine wraps the various methods for calculating upwind values.
    ! It returns a csr matrix (which must be allocated before) with the upwind
    ! values for each node pair.

    ! bucket full of fields
    type(state_type), dimension(:), intent(inout) :: state
    ! the coordinates on a mesh similar (i.e. not the same when periodic) to the field
    type(vector_field), intent(inout) :: x_field
    ! the field and its previous timelevel that we're interested in
    type(scalar_field), intent(inout) :: field, old_field
    ! the matrices of upwind and old upwind values
    type(csr_matrix), intent(inout) :: upwind_values, old_upwind_values
    ! do we want to temporarily insert the upwind elements and quadrature matrices
    ! into state so that other calls from this subroutine can use them?
    logical, optional :: defer_deletion
    ! for back compatibility pass in an option_path in case field is
    ! locally wrapped or allocated
    character(len=*), optional :: option_path

    ! a matrix containing the elements where the upwind values are projected from
    type(csr_matrix), pointer :: upwind_elements
    ! a matrix containing the quadrature within the elements where the 
    ! upwind values are projected from
    type(block_csr_matrix), pointer :: upwind_quadrature
    ! we will need the coordinates if we have to calculate upwind_elements
    type(vector_field), pointer :: x

    ! logicals controlling the type of upwind value and whether we save them
    logical :: project_point, project_grad, local, structured, l_defer_deletion, reflect, bound
    ! success indicator
    integer :: stat
    ! a local option path for back compatibility
    character(len=OPTION_PATH_LEN) :: l_option_path, spatial_discretisation_path, upwind_value_path
    ! prefix upwind matrices that have been reflected off domain boundaries
    character(len=FIELD_NAME_LEN) :: matrix_prefix
    ! which node are we projecting from?
    integer :: projection_node
    
    ewrite(2,*) 'in find_upwind_values for field ', trim(field%name), ' on sparsity from mesh ', trim(x_field%mesh%name)

    projection_node=0 ! initialise

    ! get the local option path
    if(present(option_path)) then
      l_option_path = option_path
    else
      l_option_path = trim(field%option_path)
    end if

    spatial_discretisation_path = trim(complete_cv_field_path(l_option_path))

    if(have_option(trim(spatial_discretisation_path)//"/face_value[0]/limit_face_value")) then
      upwind_value_path = trim(spatial_discretisation_path)//"/face_value[0]/limit_face_value/limiter[0]"
    else
      upwind_value_path = trim(spatial_discretisation_path)//"/face_value[0]"
    end if

    ! do we want to project to the upwind value from a point?
    project_point = have_option(trim(upwind_value_path)//&
                    '/project_upwind_value_from_point')

    ! do we want to project to the upwind value using the gradient?
    project_grad = have_option(trim(upwind_value_path)//&
                    '/project_upwind_value_from_gradient')

    ! do we want to use local values as the upwind value?
    local = have_option(trim(upwind_value_path)//&
                    '/locally_bound_upwind_value')

    ! do we want to use pseudo-structured values as the upwind value?
    structured = have_option(trim(upwind_value_path)//&
                    '/pseudo_structured_upwind_value')

    ! do we want to reflect the upwind values off the domain boundaries?
    reflect = ((have_option(trim(upwind_value_path)//&
                    '/project_upwind_value_from_point&
                    &/reflect_off_domain_boundaries')).or.&
               (have_option(trim(upwind_value_path)//&
                    '/project_upwind_value_from_gradient&
                    &/reflect_off_domain_boundaries')))
    ! do we want to bound the projected upwind values off those surrounding the upwind
    ! element?
    bound = ((have_option(trim(upwind_value_path)//&
                    '/project_upwind_value_from_point&
                    &/bound_projected_value_locally')).or.&
             (have_option(trim(upwind_value_path)//&
                    '/project_upwind_value_from_gradient&
                    &/bound_projected_value_locally')))

    ! in case none (or both) selected default to family type selection
    select case(field%mesh%shape%numbering%family)
    case (FAMILY_SIMPLEX) ! use projection
      if((.not.project_point).and.(.not.local).and.(.not.project_grad).and.(.not.structured)) then
        ewrite(2,*) "using simplex elements but haven't selected an upwind value method"
        ewrite(2,*) 'defaulting to projection from a point'
        project_point = .true.
      end if
    case (FAMILY_CUBE) ! use local
      if(project_point) then
        FLAbort("Not possible to project from a point on cube meshes")
      end if
      if(project_grad.and.bound) then
        FLAbort("Not possible to bound locally on cube meshes")
      end if
      if((.not.project_point).and.(.not.local).and.(.not.project_grad).and.(.not.structured)) then
        ewrite(2,*) "using cube elements but haven't selected an upwind value method"
        ewrite(2,*) "defaulting to a locally bound value"
        local=.true.
      end if
    case default
      FLAbort('Illegal element family')
    end select

    ! do we want to defer the deletion of the matrices
    if(present(defer_deletion)) then
      l_defer_deletion=defer_deletion
    else
      l_defer_deletion=.false.
    end if

    ! prefix the matrix name if reflected
    matrix_prefix = ""
    if(reflect) then
      matrix_prefix="Reflected"
    end if

    select case(field%field_type)
    case (FIELD_TYPE_CONSTANT)

      ! constant fields are really easy... just set all the upwind values to the constant value
      call calculate_upwind_values_constant(field, upwind_values, &
                                              old_field, old_upwind_values)

    case default
      ! do we want to project from a point?
      if(project_point) then

        upwind_elements=>extract_csr_matrix(state, &
                          trim(matrix_prefix)//trim(field%mesh%name)//"CVUpwindElements", stat)
        if(stat==0) then
        ! element matrix in state

          upwind_quadrature=>extract_block_csr_matrix(state, &
                          trim(matrix_prefix)//trim(field%mesh%name)//"CVUpwindQuadrature", stat)
          if(stat==0) then
          ! both matrices in state

              call calculate_upwind_values_project(upwind_elements, upwind_quadrature, x_field, &
                                          field, upwind_values, old_field, old_upwind_values, &
                                          bound)

          else
          ! no quadrature matrix in state

              x=>extract_vector_field(state(1), "Coordinate")

              if (l_defer_deletion.or.&
                  have_option(trim(upwind_value_path)//"/project_upwind_value_from_point&
                                  &/store_upwind_elements/store_upwind_quadrature")) then
              ! we want to save the quadrature matrix


                allocate(upwind_quadrature)
                call allocate(upwind_quadrature, upwind_values%sparsity, &
                            (/1, field%mesh%shape%loc/), &
                            name=trim(matrix_prefix)//trim(field%mesh%name)//"CVUpwindQuadrature")

                call calculate_upwind_quadrature_project(x,x_field,upwind_elements, &
                                        upwind_quadrature=upwind_quadrature, &
                                        field=field, upwind_values=upwind_values, &
                                        old_field=old_field, old_upwind_values=old_upwind_values, &
                                        reflect=reflect, bound=bound)

                call insert(state, upwind_quadrature, trim(upwind_quadrature%name))

                call deallocate(upwind_quadrature)
                deallocate(upwind_quadrature)

              else
              ! don't even calculate the quadrature matrix

                call calculate_upwind_quadrature_project(x,x_field,upwind_elements, &
                                        field=field, upwind_values=upwind_values, &
                                        old_field=old_field, old_upwind_values=old_upwind_values, &
                                        reflect=reflect, bound=bound)

              end if
          end if
        else
        ! no element matrix in state, may as well make both now

          x=>extract_vector_field(state(1), "Coordinate")

          if (l_defer_deletion.or.&
              have_option(trim(upwind_value_path)//"/project_upwind_value_from_point&
                                &/store_upwind_elements")) then
          ! we want to save the element matrix

              allocate(upwind_elements)
              call allocate(upwind_elements, upwind_values%sparsity, &
                          type=CSR_INTEGER, &
                          name=trim(matrix_prefix)//trim(field%mesh%name)//"CVUpwindElements")

              if (l_defer_deletion.or.&
                  have_option(trim(upwind_value_path)//"/project_upwind_value_from_point&
                                  &/store_upwind_elements/store_upwind_quadrature")) then
              ! we want to save both matrices

                allocate(upwind_quadrature)
                call allocate(upwind_quadrature, upwind_values%sparsity, &
                            (/1, field%mesh%shape%loc/), &
                            name=trim(matrix_prefix)//trim(field%mesh%name)//"CVUpwindQuadrature")

                call calculate_all_upwind_project(x, x_field, upwind_elements=upwind_elements, &
                                        upwind_quadrature=upwind_quadrature, &
                                        field=field, upwind_values=upwind_values, &
                                        old_field=old_field, old_upwind_values=old_upwind_values, &
                                        reflect=reflect, bound=bound)

                call insert(state, upwind_quadrature, trim(upwind_quadrature%name))

                call deallocate(upwind_quadrature)
                deallocate(upwind_quadrature)

              else
              ! we want to save the element matrix but not the quadrature

                call calculate_all_upwind_project(x, x_field, upwind_elements=upwind_elements, &
                                        field=field, upwind_values=upwind_values, &
                                        old_field=old_field, old_upwind_values=old_upwind_values, &
                                        reflect=reflect, bound=bound)

              end if

              call insert(state, upwind_elements, trim(upwind_elements%name))

              call deallocate(upwind_elements)
              deallocate(upwind_elements)

          else
          ! we don't want anything but the values

              call calculate_all_upwind_project(x, x_field, field=field, upwind_values=upwind_values, &
                                        old_field=old_field, old_upwind_values=old_upwind_values, &
                                        reflect=reflect, bound=bound)

          end if

        end if

      ! we don't want to project from a point, shall we use the gradient instead?
      else if(project_grad) then

        projection_node=cv_projection_node(trim(l_option_path))
        x=>extract_vector_field(state(1), "Coordinate")

        if(bound) then
          ! we need to know about upwind_elements so we can bound... do we have them in state?
          upwind_elements=>extract_csr_matrix(state, &
                            trim(matrix_prefix)//trim(field%mesh%name)//"CVUpwindElements", stat)
          if(stat==0) then
          ! element matrix in state
             call calculate_upwind_values_project_grad(x, x_field, upwind_elements=upwind_elements, &
                                                       field=field, upwind_values=upwind_values, &
                                                       old_field=old_field, old_upwind_values=old_upwind_values, &
                                                       reflect=reflect, bound=bound, projection_node=projection_node)

          else
            ! no element matrix in state... do we want to save it for next time?
            if (l_defer_deletion.or.&
               have_option(trim(upwind_value_path)//"/project_upwind_value_from_gradient&
                                 &/bound_projected_value_locally/store_upwind_elements")) then
               ! yes
               allocate(upwind_elements)
               call allocate(upwind_elements, upwind_values%sparsity, &
                           type=CSR_INTEGER, &
                           name=trim(matrix_prefix)//trim(field%mesh%name)//"CVUpwindElements")

               call calculate_all_upwind_project_grad(x, x_field, upwind_elements=upwind_elements, &
                                                      field=field, upwind_values=upwind_values, &
                                                      old_field=old_field, old_upwind_values=old_upwind_values, &
                                                      reflect=reflect, bound=bound, projection_node=projection_node)

               call insert(state, upwind_elements, trim(upwind_elements%name))

               call deallocate(upwind_elements)
               deallocate(upwind_elements)

            else
               ! no
               call calculate_all_upwind_project_grad(x, x_field, &
                                                      field=field, upwind_values=upwind_values, &
                                                      old_field=old_field, old_upwind_values=old_upwind_values, &
                                                      reflect=reflect, bound=bound, projection_node=projection_node)

            end if
          end if

        else

          ! we're not bounding so there's no need to worry about upwind_elements
          call calculate_all_upwind_project_grad(x, x_field, field=field, upwind_values=upwind_values, &
                                                 old_field=old_field, old_upwind_values=old_upwind_values, &
                                                 reflect=reflect, bound=bound, projection_node=projection_node)

        end if

      ! ok, we don't want to project at all, so do we want to use local values?
      else if(local) then

        call calculate_upwind_values_local(field, upwind_values, old_field, old_upwind_values)

      ! a completely new way... let's try pseudo-structured
      else if(structured) then

        call calculate_upwind_values_structured(x_field, field, upwind_values, &
                                           old_field, old_upwind_values)

      ! no, um, something's gone wrong here then!
      else

        FLAbort("Unknown upwind value calculation method.")

      end if

    end select

  end subroutine find_upwind_values_multiple_states

  subroutine calculate_upwind_values_constant(field, upwind_values, &
                                              old_field, old_upwind_values)

      ! just set the upwind values to the constant field values
      ! should only be called for FIELD_TYPE_CONSTANT fields

      type(scalar_field), intent(in) :: field, old_field
      type(csr_matrix), intent(inout) :: upwind_values, old_upwind_values

      upwind_values%val = field%val(1)
      old_upwind_values%val = old_field%val(1)

   end subroutine calculate_upwind_values_constant

  subroutine calculate_all_upwind_project(x, x_field, upwind_elements, upwind_quadrature, &
                                  field, upwind_values, old_field, old_upwind_values, &
                                  reflect, bound)

    ! project from a node pair to an upwind value when we have no information available at
    ! all... i.e. we need to calculate which element the upwind value is in, then we need to 
    ! work out its quadrature, then we need to actually find the value and possibly bound it.

    ! coordinates
    type(vector_field), intent(inout) :: x, x_field

    type(csr_matrix), intent(inout), optional :: upwind_elements
    type(block_csr_matrix), intent(inout), optional :: upwind_quadrature

    type(scalar_field), intent(inout) :: field
    type(scalar_field), intent(in) :: old_field
    type(csr_matrix), intent(inout) :: upwind_values
    type(csr_matrix), intent(inout) :: old_upwind_values
    logical, intent(in) :: reflect, bound

    ! local memory:
    ! allocatable memory
    integer, dimension(:), pointer :: nodes, eles, x_eles
    ! loop integers
    integer :: i, j, k, ele, l_ele, dim, local_coord, i_field

    ! the coordinates of a pt just upwind of the node pair
    real, dimension(mesh_dim(field)) :: xc, xc_vector
    ! element node coordinates
    real, dimension(mesh_dim(x), ele_loc(x, 1)) :: x_ele
    real, dimension(mesh_dim(x_field), ele_loc(x_field, 1)) :: x_field_ele
    ! simplex volume coordinates
    real, dimension(field%mesh%shape%quadrature%vertices) :: coords, l_coords
    ! element field node values
    real, dimension(field%mesh%shape%loc) :: field_ele
    ! the global node numbers of an element
    integer, dimension(:), pointer :: field_nodes
    ! the upwind value
    real :: upwind_value

    real, dimension(field%mesh%shape%loc) :: l_shape

    ! a vector field of normals if we're reflecting
    type(vector_field) :: normals
    ! a logical list saying which nodes are on the boundary
    logical, dimension(:), allocatable :: on_boundary
    ! which bc are on which nodes
    integer, dimension(:), allocatable :: field_bc_type

    logical :: upwind_elements_present, upwind_quadrature_present

    integer :: save_pos=0 ! saves the position in the matrix for optimisation
    
    ewrite(2,*) 'in calculate_all_upwind_project'
    ! the projected point values upwind value matrix is on the x_field mesh
    ! which cannot be periodic

    upwind_elements_present=.false.
    upwind_quadrature_present=.false.

    ! zero everything we have
    call zero(upwind_values)
    call zero(old_upwind_values)
    if(present(upwind_elements)) then
      call zero(upwind_elements)
      upwind_elements_present=.true.
    end if
    if(present(upwind_quadrature)) then
      call zero(upwind_quadrature)
      upwind_quadrature_present=.true.
    end if

    call allocate(normals, mesh_dim(x_field), x_field%mesh, name="NormalsToBoundary")
    call zero(normals)
    allocate(on_boundary(node_count(x_field)))
    on_boundary=.false.
    if(reflect) then
      ! work out what the domain normals are
      call calculate_boundary_normals(field%mesh, x, &
                                      normals, on_boundary)
    end if

    allocate(field_bc_type(node_count(field)))
    field_bc_type = 0
    ! create the node to element list
    call add_nelist(x_field%mesh)
    if(mesh_periodic(field)) then
      call add_nelist(field%mesh)
      call get_boundary_condition_nodes(field, (/"periodic"/), field_bc_type)
    end if

    dim = mesh_dim(field)
    if((dim/=2).and.(dim/=3)) then
      FLAbort("Unsupported dimension count")
    end if

    ! loop over the nodes
    do i = 1, size(upwind_values,1)
      ! find the neighbouring nodes using the matrix sparsity (of the unperiodic coordinate mesh!)
      ! FIXME: the matrix sparsity is not the same as the mesh connectivity for 
      !        cube meshes!
      !        (not very important for projection as it only works on simplex meshes)
      nodes => row_m_ptr(upwind_values, i)
      if (size(nodes) == 0) cycle
      ! find the neighbouring elements using the node to element list
      x_eles => node_neigh(x_field, i)
      if(mesh_periodic(field)) then
        local_coord = local_coords(x_field, x_eles(1), i)
        field_nodes=>ele_nodes(field, x_eles(1))
        i_field = field_nodes(local_coord)
        eles => node_neigh(field, i_field)
      else
        i_field = i
        eles => x_eles
      end if
      ! loop over the neighbouring nodes
      do j = 1, size(nodes)
        if(nodes(j)==i) cycle  ! skip the node that's the same as the i node

        ! find the vector connecting to the point just upwind of the node pair
        ! (also deals with reflection)
        xc_vector=project_upwind(i, nodes(j), x_field, &
                          on_boundary, normals)
                          
        if(field_bc_type(i_field)==0) then
          xc = node_val(x_field, i) + xc_vector
        end if

        l_coords = infinity
        l_ele = eles(1)

        ! loop over neighbouring elements working out which one contains
        ! (or nearly contains) xc
        do k = 1, size(eles)
          ele = eles(k)
          x_ele=ele_val(x, ele)
          
          if(field_bc_type(i_field)==1) then
            ! find the local node number so we can add the vector
            ! pointing at the upwind point to the coordinates at
            ! that node
            ! (this is done this way in case we're periodic and on
            ! a boundary)
            x_field_ele=ele_val(x_field, ele)
            local_coord = local_coords(field, ele, i_field)
            xc = x_field_ele(:,local_coord) + xc_vector
          end if

          select case(dim)
          case(2)
            coords=calculate_area_coordinates(x_ele, xc) ! only makes sense with linear coordinates in x_ele
          case(3)
            coords=calculate_volume_coordinates(x_ele, xc) ! only makes sense with linear coordinates in x_ele
          end select

          if(sum(coords)<sum(l_coords)) then
            l_coords=coords  ! save the one we think contains xc
            l_ele=ele        !
          end if

        end do
        
        ! just in case (for instance when we haven't reflected off a domain boundary)
        ! make sure all our quadratures are positive
        do k = 1, size(l_coords)
          l_coords(k) = max(0.0, l_coords(k))
        end do
        ! normalise
        l_coords = l_coords/sum(l_coords)

        l_shape = eval_shape(field%mesh%shape, l_coords)

        ! if we want them then assemble the upwind_elements and upwind_quadrature matrices
        if(upwind_elements_present) then
          call set(upwind_elements, i, nodes(j), l_ele, save_pos=save_pos)
        end if
        if(upwind_quadrature_present) then
           do k = 1, size(l_shape)
              call set(upwind_quadrature, 1, k, i, nodes(j), l_shape(k), save_pos=save_pos)
           end do
        end if

        ! calculate the upwind value:
        field_ele=ele_val(field, l_ele)
        upwind_value=dot_product(l_shape,field_ele)
        ! project it a little further
        upwind_value=node_val(field, i_field) + (1./c_distance)*(upwind_value-node_val(field, i_field))
        if(bound) then
           ! bound it relative to the surrounding values
           upwind_value = max(min(upwind_value, maxval(field_ele)), minval(field_ele))
        end if
        ! set in the matrix
        call set(upwind_values, i, nodes(j), upwind_value, save_pos=save_pos)

        ! same for the old field
        field_ele=ele_val(old_field, l_ele)
        upwind_value=dot_product(l_shape,field_ele)
        upwind_value=node_val(old_field, i_field) + (1./c_distance)*(upwind_value-node_val(old_field, i_field))
        if(bound) then
            upwind_value = max(min(upwind_value, maxval(field_ele)), minval(field_ele))
        end if
        call set(old_upwind_values, i, nodes(j), upwind_value, save_pos=save_pos)

      end do
    end do

    call deallocate(normals)
    deallocate(field_bc_type)

  end subroutine calculate_all_upwind_project

  subroutine calculate_upwind_quadrature_project(x,x_field,upwind_elements, upwind_quadrature, &
                  field, upwind_values, old_field, old_upwind_values, &
                  reflect, bound)
    ! project from a node pair to an upwind value when we have the upwind_elements available 
    ! ... i.e. we need to work out each elements its quadrature, then we need to actually 
    ! find the value and possibly bound it.

    ! coordinates
    type(vector_field), intent(inout) :: x, x_field
    type(csr_matrix), intent(in) :: upwind_elements
    type(scalar_field), intent(inout) :: field
    type(scalar_field), intent(in) :: old_field
    type(csr_matrix), intent(inout) :: upwind_values
    type(csr_matrix), intent(inout) :: old_upwind_values

    type(block_csr_matrix), intent(inout), optional :: upwind_quadrature
    logical, intent(in) :: reflect, bound

    ! local memory:
    ! integer loops
    integer :: i, j, k, l_ele, dim, local_coord, i_field
    ! allocatable memory
    integer, dimension(:), pointer :: nodes, x_eles

    ! the coordinates of a pt just upwind of the node pair
    real, dimension(x%dim) :: xc, xc_vector
    real, dimension(x%dim, ele_loc(x, 1)) :: x_ele
    real, dimension(x_field%dim, ele_loc(x_field, 1)) :: x_field_ele
    real, dimension(field%mesh%shape%quadrature%vertices) :: l_coords
    real, dimension(field%mesh%shape%loc) :: field_ele
    integer, dimension(:), pointer :: field_nodes
    real :: upwind_value

    real, dimension(field%mesh%shape%loc) :: l_shape

    type(vector_field) :: normals
    logical, dimension(:), allocatable :: on_boundary
    ! which bc are on which nodes
    integer, dimension(:), allocatable :: field_bc_type

    logical :: upwind_quadrature_present

    integer :: save_pos=0 ! saves the position in the matrix for optimisation
    
    ewrite(2,*) 'in calculate_upwind_quadrature_project'
    ! the projected point values upwind value matrix is on the x_field mesh
    ! which cannot be periodic

    upwind_quadrature_present=.false.

    ! zero everything we have
    call zero(upwind_values)
    call zero(old_upwind_values)
    if(present(upwind_quadrature)) then
      call zero(upwind_quadrature)
      upwind_quadrature_present=.true.
    end if

    dim = mesh_dim(field)
    if((dim/=2).and.(dim/=3)) then
      FLAbort("Unsupported dimension count")
    end if

    call allocate(normals, mesh_dim(x_field), x_field%mesh, name="NormalsToBoundary")
    call zero(normals)
    allocate(on_boundary(node_count(x_field)))
    on_boundary=.false.
    if(reflect) then
      call calculate_boundary_normals(field%mesh, x, &
                                      normals, on_boundary)
    end if
    
    allocate(field_bc_type(node_count(field)))
    field_bc_type = 0
    if(mesh_periodic(field)) then
      call add_nelist(x_field%mesh)
      call get_boundary_condition_nodes(field, (/"periodic"/), field_bc_type)
    end if

    do i = 1, size(upwind_values, 1)
      ! find the neighbouring nodes using the matrix sparsity of the unperiodic mesh
      ! FIXME: the matrix sparsity is not the same as the mesh connectivity for 
      !        cube meshes!
      !        (not very important for projection as it only works on simplex meshes)
        nodes => row_m_ptr(upwind_values, i)
        
        if(mesh_periodic(field)) then
          x_eles=>node_neigh(x_field, i)
          local_coord = local_coords(x_field, x_eles(1), i)
          field_nodes=>ele_nodes(field, x_eles(1))
          i_field = field_nodes(local_coord)
        else
          i_field = i
        end if
        do j = 1, size(nodes)
          if(nodes(j)==i) cycle
          
          ! find the vector connecting to the point just upwind of the node pair
          ! (also deals with reflection)
          xc_vector=project_upwind(i, nodes(j), x_field, &
                                   on_boundary, normals)

          l_ele=ival(upwind_elements,i,nodes(j))
          x_ele=ele_val(x, l_ele)
          
          if(field_bc_type(i_field)==1) then
            ! find the local node number so we can add the vector
            ! pointing at the upwind point to the coordinates at
            ! that node
            ! (this is done this way in case we're periodic and on
            ! a boundary)
            x_field_ele=ele_val(x_field, l_ele)
            local_coord = local_coords(field, l_ele, i_field)
            xc = x_field_ele(:,local_coord) + xc_vector
          else
            xc = node_val(x_field, i) + xc_vector
          end if
          
          select case(dim)
          case(2)
              l_coords=calculate_area_coordinates(x_ele, xc) ! only makes sense with linear coordinates in x_ele
          case(3)
              l_coords=calculate_volume_coordinates(x_ele, xc) ! only makes sense with linear coordinates in x_ele
          end select

          do k = 1, size(l_coords)
              l_coords(k) = max(0.0, l_coords(k))
          end do
          l_coords = l_coords/sum(l_coords)

          l_shape = eval_shape(field%mesh%shape, l_coords)

          if(upwind_quadrature_present) then
              do k = 1, size(l_shape)
                call set(upwind_quadrature, 1, k, i, nodes(j), l_shape(k), save_pos=save_pos)
              end do
          end if

          field_ele=ele_val(field, l_ele)
          upwind_value=dot_product(l_shape,field_ele)
          upwind_value=node_val(field, i_field) + (1./c_distance)*(upwind_value-node_val(field, i_field))
          if(bound) then
              upwind_value = max(min(upwind_value, maxval(field_ele)), minval(field_ele))
          end if
          call set(upwind_values, i, nodes(j), upwind_value, save_pos=save_pos)

          field_ele=ele_val(old_field, l_ele)
          upwind_value=dot_product(l_shape,field_ele)
          upwind_value=node_val(old_field, i_field) + (1./c_distance)*(upwind_value-node_val(old_field, i_field))
          if(bound) then
              upwind_value = max(min(upwind_value, maxval(field_ele)), minval(field_ele))
          end if
          call set(old_upwind_values, i, nodes(j), upwind_value, save_pos=save_pos)

        end do
    end do

    call deallocate(normals)
    deallocate(field_bc_type)

  end subroutine calculate_upwind_quadrature_project

  subroutine calculate_upwind_values_project(upwind_elements, upwind_quadrature, x_field, &
                                             field, upwind_values, &
                                             old_field, old_upwind_values, &
                                             bound)
      ! project from a node pair to an upwind value when we have the upwind_elements and
      ! the upwind_quadrature available... i.e. we only need to
      ! find the value and possibly bound it.

      type(csr_matrix), intent(in) :: upwind_elements
      type(block_csr_matrix), intent(in) :: upwind_quadrature
      
      type(vector_field), intent(inout) :: x_field

      type(scalar_field), intent(in) :: field
      type(scalar_field), intent(in) :: old_field
      type(csr_matrix), intent(inout) :: upwind_values
      type(csr_matrix), intent(inout) :: old_upwind_values

      ! the global node numbers of an element
      integer, dimension(:), pointer :: field_nodes, x_eles

      logical, intent(in) :: bound

      integer :: i, j, k, l_ele, local_coord, i_field
      integer, dimension(:), pointer :: nodes
      real, dimension(field%mesh%shape%loc) :: field_ele
      real, dimension(field%mesh%shape%loc) :: l_shape
      real :: upwind_value

      integer :: save_pos=0 ! saves the position in the matrix for optimisation

      ewrite(2,*) 'in calculate_upwind_values_project'
      ! the projected point values upwind value matrix is on the x_field mesh
      ! which cannot be periodic

      call zero(upwind_values)
      call zero(old_upwind_values)

      if(mesh_periodic(field)) then
        call add_nelist(x_field%mesh)
      end if

      do i = 1, size(upwind_values, 1)
        ! find the neighbouring nodes using the matrix sparsity
        ! FIXME: the matrix sparsity is not the same as the mesh connectivity for 
        !        cube meshes!
        !        (not very important for projection as it only works on simplex meshes)
         nodes => row_m_ptr(upwind_values, i)
       
         if (size(nodes) == 0) cycle

         if(mesh_periodic(field)) then
           x_eles=>node_neigh(x_field, i)
           local_coord = local_coords(x_field, x_eles(1), i)
           field_nodes=>ele_nodes(field, x_eles(1))
           i_field = field_nodes(local_coord)
         else
           i_field = i
         end if
         do j = 1, size(nodes)
            if(nodes(j)==i) cycle

            do k = 1, size(l_shape)
               l_shape(k)=val(upwind_quadrature, 1, k, i, nodes(j), save_pos=save_pos)
            end do

            l_ele = ival(upwind_elements, i, nodes(j))
          
            field_ele=ele_val(field, l_ele)
            upwind_value=dot_product(l_shape,field_ele)
            upwind_value=node_val(field, i_field) + (1./c_distance)*(upwind_value-node_val(field, i_field))
            if(bound) then
               upwind_value = max(min(upwind_value, maxval(field_ele)), minval(field_ele))
            end if
            call set(upwind_values, i, nodes(j), upwind_value, save_pos=save_pos)

            field_ele=ele_val(old_field, ival(upwind_elements, i, nodes(j)))
            upwind_value=dot_product(l_shape,field_ele)
            upwind_value=node_val(old_field, i_field) + (1./c_distance)*(upwind_value-node_val(old_field, i_field))
            if(bound) then
                upwind_value = max(min(upwind_value, maxval(field_ele)), minval(field_ele))
            end if
            call set(old_upwind_values, i, nodes(j), upwind_value, save_pos=save_pos)

         end do
      end do

   end subroutine calculate_upwind_values_project

  subroutine calculate_all_upwind_project_grad(x, x_field, upwind_elements, &
                                  field, upwind_values, old_field, old_upwind_values, &
                                  reflect, bound, projection_node)

    ! project from a node pair to an upwind value using the interpolated gradient of the fi.

    ! coordinates
    type(vector_field), intent(inout) :: x, x_field

    type(csr_matrix), intent(inout), optional :: upwind_elements

    type(scalar_field), intent(inout) :: field
    type(scalar_field), intent(in) :: old_field
    type(csr_matrix), intent(inout) :: upwind_values
    type(csr_matrix), intent(inout) :: old_upwind_values
    logical, intent(in) :: reflect, bound
    integer, intent(in) :: projection_node

    ! local memory:
    ! allocatable memory
    integer, dimension(:), pointer :: nodes, eles, x_eles, field_nodes
    ! loop integers
    integer :: i, j, k, ele, l_ele, dim, local_coord, i_field, j_field

    ! the coordinates of a pt just upwind of the node pair and the vector
    ! in the direction of projection
    real, dimension(x%dim) :: xc, xc_vector, d
    ! the gradients of the field and old_field at the donor node
    real, dimension(x%dim) :: grad_c, old_grad_c
    ! element node coordinates
    real, dimension(x%dim, ele_loc(x, 1)) :: x_ele
    real, dimension(x_field%dim, ele_loc(x_field, 1)) :: x_field_ele
    ! simplex volume coordinates
    real, dimension(x%mesh%shape%loc) :: coords, l_coords
    ! element field node values
    real, dimension(field%mesh%shape%loc) :: field_ele
    ! the upwind value
    real :: upwind_value, old_upwind_value

    ! a vector field of normals if we're reflecting
    type(vector_field) :: normals
    ! a logical list saying which nodes are on the boundary
    logical, dimension(:), allocatable :: on_boundary
    ! which bc are on which nodes
    integer, dimension(:), allocatable :: field_bc_type

    ! gradients
    type(vector_field) :: grad_field, grad_old_field

    logical :: upwind_elements_present

    integer :: save_pos=0 ! saves the position in the matrix for optimisation

    ewrite(2,*) 'in calculate_all_upwind_project_grad'
    ! the projected gradient values upwind value matrix is on the x_field mesh
    ! which cannot be periodic
    
    upwind_elements_present=.false.

    ! zero everything we have
    call zero(upwind_values)
    call zero(old_upwind_values)
    if(present(upwind_elements)) then
      call zero(upwind_elements)
      upwind_elements_present=.true.
    end if

    dim = mesh_dim(field)
    if((dim/=2).and.(dim/=3)) then
      FLAbort("Unsupported dimension count")
    end if

    call allocate(normals, dim, x_field%mesh, name="NormalsToBoundary")
    call zero(normals)
    allocate(on_boundary(node_count(x_field)))
    on_boundary=.false.
    if(reflect) then
      ! work out what the domain normals are
      call calculate_boundary_normals(field%mesh, x, &
                                      normals, on_boundary)
    end if

    call allocate(grad_field, dim, field%mesh, name="FieldGradient")
    call zero(grad_field)
    call grad(field, x, grad_field)

    call allocate(grad_old_field, dim, field%mesh, name="OldFieldGradient")
    call zero(grad_old_field)
    call grad(old_field, x, grad_old_field)

    allocate(field_bc_type(node_count(field)))
    field_bc_type=0
    if(bound.or.mesh_periodic(field)) then
      ! create the node to element list
      call add_nelist(x_field%mesh)
      if(bound.and.mesh_periodic(field)) then
        call add_nelist(field%mesh)
        call get_boundary_condition_nodes(field, (/"periodic"/), field_bc_type)
      end if
    end if
    eles => null()
    xc = 0.0
    xc_vector = 0.0
    
    ! loop over the nodes
    do i = 1, size(upwind_values,1)
      ! find the neighbouring nodes using the matrix sparsity
      ! FIXME: the matrix sparsity is not the same as the mesh connectivity for 
      !        cube meshes!
      nodes => row_m_ptr(upwind_values, i)
      if(bound.or.mesh_periodic(field)) then
         ! find the neighbouring elements using the node to element list
        x_eles => node_neigh(x_field, i)
        if(mesh_periodic(field)) then
          local_coord = local_coords(x_field, x_eles(1), i)
          field_nodes=>ele_nodes(field, x_eles(1))
          i_field = field_nodes(local_coord)
          if(bound) then
            eles => node_neigh(field, i_field)
          end if
        else
          i_field = i
          if(bound) then
            eles => x_eles
          end if
        end if
      else
        i_field = i
      end if
      ! loop over the neighbouring nodes
      do j = 1, size(nodes)
        if(nodes(j)==i) cycle  ! skip the node that's the same as the i node
        ! i is considered as the donor node and j is considered to be the downwind
        ! (although we don't actually know this yet)

         d=project_vector(i, nodes(j), x_field, &
                          on_boundary, normals) ! d is the vector between the donor and downwind nodes

         grad_c = node_val(grad_field, i_field)          ! the gradient at the donor node
         old_grad_c = node_val(grad_old_field, i_field)  ! the old gradient at the donor node

         if(bound) then

            ! find the vecotr connecting to the point just upwind of the node pair
            ! (also deals with reflection)
            xc_vector=project_upwind(i, nodes(j), x_field, &
                              on_boundary, normals)

            if(field_bc_type(i_field)==0) then
              xc = node_val(x_field, i) + xc_vector
            end if

            l_coords = infinity
            l_ele = eles(1)

            ! loop over neighbouring elements working out which one contains
            ! (or nearly contains) xc
            do k = 1, size(eles)
               ele = eles(k)
               x_ele=ele_val(x, ele)
          
               if(field_bc_type(i_field)==1) then
                 ! find the local node number so we can add the vector
                 ! pointing at the upwind point to the coordinates at
                 ! that node
                 ! (this is done this way in case we're periodic and on
                 ! a boundary)
                 x_field_ele=ele_val(x_field, ele)
                 local_coord = local_coords(field, ele, i_field)
                 xc = x_field_ele(:,local_coord) + xc_vector
               end if

               select case(dim)
               case(2)
                  coords=calculate_area_coordinates(x_ele, xc) ! only makes sense with linear coordinates in x_ele
               case(3)
                  coords=calculate_volume_coordinates(x_ele, xc) ! only makes sense with linear coordinates in x_ele
               end select

               if(sum(coords)<sum(l_coords)) then
                  l_coords=coords  ! save the one we think contains xc
                  l_ele=ele        !
               end if
            end do

            if(upwind_elements_present) then
               call set(upwind_elements, i, nodes(j), l_ele, save_pos=save_pos)
            end if

         end if

         select case(projection_node)
         case(CV_DOWNWIND_PROJECTION_NODE)
            if(mesh_periodic(field)) then
              x_eles => node_neigh(x_field, nodes(j))
              local_coord = local_coords(x_field, x_eles(1), nodes(j))
              field_nodes=>ele_nodes(field, x_eles(1))
              j_field = field_nodes(local_coord)
            else
              j_field = j
            end if
            ! project from downwind node (Jasak et al., 1999)
            upwind_value=node_val(field, j_field)-2.0*dot_product(d, grad_c)
            old_upwind_value=node_val(old_field, j_field)-2.0*dot_product(d, old_grad_c)
         case(CV_DONOR_PROJECTION_NODE)
            ! or project from donor node
            upwind_value=node_val(field, i_field)-dot_product(d, grad_c)
            old_upwind_value=node_val(old_field, i_field)-dot_product(d, old_grad_c)
         end select

         if(bound) then
            ! calculate neighbouring values:
            field_ele=ele_val(field, l_ele)
            ! bound it relative to the surrounding values
            upwind_value = max(min(upwind_value, maxval(field_ele)), minval(field_ele))

            ! and for the old field
            field_ele=ele_val(old_field, l_ele)
            ! bound it relative to the surrounding values
            old_upwind_value = max(min(old_upwind_value, maxval(field_ele)), minval(field_ele))
         end if
         ! set in the matrix
         call set(upwind_values, i, nodes(j), upwind_value, save_pos=save_pos)
         call set(old_upwind_values, i, nodes(j), old_upwind_value, save_pos=save_pos)

      end do
    end do

    call deallocate(normals)
    call deallocate(grad_field)
    call deallocate(grad_old_field)
    deallocate(field_bc_type)

  end subroutine calculate_all_upwind_project_grad

  subroutine calculate_upwind_values_project_grad(x, x_field, upwind_elements, &
                                  field, upwind_values, old_field, old_upwind_values, &
                                  reflect, bound, projection_node)

    ! project from a node pair to an upwind value using the interpolated gradient of the fi.

    ! coordinates
    type(vector_field), intent(inout) :: x, x_field

    type(csr_matrix), intent(in) :: upwind_elements

    type(scalar_field), intent(inout) :: field
    type(scalar_field), intent(in) :: old_field
    type(csr_matrix), intent(inout) :: upwind_values
    type(csr_matrix), intent(inout) :: old_upwind_values
    logical, intent(in) :: reflect, bound
    integer, intent(in) :: projection_node

    ! local memory:
    ! allocatable memory
    integer, dimension(:), pointer :: nodes, field_nodes, x_eles
    ! loop integers
    integer :: i, j, l_ele, dim, local_coord, i_field, j_field

    ! the vector in the direction of projection
    real, dimension(x%dim) :: d
    ! the gradients of the field and old_field at the donor node
    real, dimension(x%dim) :: grad_c, old_grad_c
    ! element field node values
    real, dimension(field%mesh%shape%loc) :: field_ele
    ! the upwind value
    real :: upwind_value, old_upwind_value

    ! a vector field of normals if we're reflecting
    type(vector_field) :: normals
    ! a logical list saying which nodes are on the boundary
    logical, dimension(:), allocatable :: on_boundary

    ! gradients
    type(vector_field) :: grad_field, grad_old_field

    integer :: save_pos=0 ! saves the position in the matrix for optimisation

    ewrite(2,*) 'in calculate_upwind_values_project_grad'
    ! the projected gradient values upwind value matrix is on the x_field mesh
    ! which cannot be periodic

    ! zero everything we have
    call zero(upwind_values)
    call zero(old_upwind_values)

    dim = mesh_dim(field)
    if((dim/=2).and.(dim/=3)) then
      FLAbort("Unsupported dimension count")
    end if

    call allocate(normals, dim, x_field%mesh, name="NormalsToBoundary")
    call zero(normals)
    allocate(on_boundary(node_count(x_field)))
    on_boundary=.false.
    if(reflect) then
      ! work out what the domain normals are
      call calculate_boundary_normals(field%mesh, x, &
                                      normals, on_boundary)
    end if

    if(mesh_periodic(field)) then
      call add_nelist(x_field%mesh)
    end if

    call allocate(grad_field, dim, field%mesh, name="FieldGradient")
    call zero(grad_field)
    call grad(field, x, grad_field)

    call allocate(grad_old_field, dim, field%mesh, name="OldFieldGradient")
    call zero(grad_old_field)
    call grad(old_field, x, grad_old_field)

    ! loop over the nodes
    do i = 1, size(upwind_values,1)
      ! find the neighbouring nodes using the matrix sparsity
      ! FIXME: the matrix sparsity is not the same as the mesh connectivity for 
      !        cube meshes!
      !        (not very important for projection as it only works on simplex meshes)
      nodes => row_m_ptr(upwind_values, i)
      if(mesh_periodic(field)) then
         ! find the neighbouring elements using the node to element list
        x_eles => node_neigh(x_field, i)
        local_coord = local_coords(x_field, x_eles(1), i)
        field_nodes=>ele_nodes(field, x_eles(1))
        i_field = field_nodes(local_coord)
      else
        i_field = i
      end if

      ! loop over the neighbouring nodes
      do j = 1, size(nodes)
        if(nodes(j)==i) cycle  ! skip the node that's the same as the i node
        ! i is considered as the donor node and j is considered to be the downwind
        ! (although we don't actually know this yet)

         d=project_vector(i, nodes(j), x_field, &
                          on_boundary, normals) ! d is the vector between the donor and downwind nodes

         grad_c = node_val(grad_field, i_field)          ! the gradient at the donor node
         old_grad_c = node_val(grad_old_field, i_field)  ! the old gradient at the donor node

         select case(projection_node)
         case(CV_DOWNWIND_PROJECTION_NODE)
            ! project from downwind node (Jasak et al., 1999)
            if(mesh_periodic(field)) then
              x_eles => node_neigh(x_field, nodes(j))
              local_coord = local_coords(x_field, x_eles(1), nodes(j))
              field_nodes=>ele_nodes(field, x_eles(1))
              j_field = field_nodes(local_coord)
            else
              j_field = j
            end if
            upwind_value=node_val(field, j_field)-2.0*dot_product(d, grad_c)
            old_upwind_value=node_val(old_field, j_field)-2.0*dot_product(d, old_grad_c)
         case(CV_DONOR_PROJECTION_NODE)
            ! or project from donor node
            upwind_value=node_val(field, i_field)-dot_product(d, grad_c)
            old_upwind_value=node_val(old_field, i_field)-dot_product(d, old_grad_c)
         end select

         if(bound) then
            ! which element are we in?
            l_ele=ival(upwind_elements,i,nodes(j))

            ! calculate neighbouring values:
            field_ele=ele_val(field, l_ele)
            ! bound it relative to the surrounding values
            upwind_value = max(min(upwind_value, maxval(field_ele)), minval(field_ele))

            ! and for the old field
            field_ele=ele_val(old_field, l_ele)
            ! bound it relative to the surrounding values
            old_upwind_value = max(min(old_upwind_value, maxval(field_ele)), minval(field_ele))
         end if
         ! set in the matrix
         call set(upwind_values, i, nodes(j), upwind_value, save_pos=save_pos)
         call set(old_upwind_values, i, nodes(j), old_upwind_value, save_pos=save_pos)

      end do
    end do

    call deallocate(normals)
    call deallocate(grad_field)
    call deallocate(grad_old_field)
    
  end subroutine calculate_upwind_values_project_grad

  subroutine calculate_upwind_values_local(field, upwind_values, &
                                           old_field, old_upwind_values)
      ! use the local values surrounding a node pair to calculate the
      ! upwind value (choose the min or max depending on the gradient
      ! between the node pair)

      type(scalar_field), intent(in) :: field
      type(scalar_field), intent(in) :: old_field
      type(csr_matrix), intent(inout) :: upwind_values
      type(csr_matrix), intent(inout) :: old_upwind_values

      integer :: i, j
      integer, dimension(:), pointer :: nodes

      integer :: save_pos=0 ! saves the position in the matrix for optimisation
                 ! this is of dubious benefit here as the value setting is potentially
                 ! different for current and old values
      ewrite(2,*) 'in calculate_upwind_values_local'
      ! the local values upwind value matrix is on the field mesh
      ! which may be periodic
      ! (this is done because it doesn't depend on coordinates like the other methods)

      ! zero everything we have
      call zero(upwind_values)
      call zero(old_upwind_values)
      
      ! loop over nodes
      do i = 1, size(upwind_values, 1)
        ! find the neighbouring nodes using the matrix sparsity
        ! FIXME: the matrix sparsity is not the same as the mesh connectivity for 
        !        cube meshes!
        !        (this could potentially be very important for local upwind values
        !         as unnconected values may be used)
         nodes => row_m_ptr(upwind_values, i)
         ! loop over neighbouring nodes
         do j = 1, size(nodes)

            ! don't try using the donor as the downwind
            if(nodes(j)==i) cycle

            ! test gradient between nodes and calculate min or max of neighbouring node values
            if(node_val(field, i)<node_val(field, nodes(j))) then
              call set(upwind_values, i, nodes(j), minval(field%val(nodes), mask=((nodes/=i).and.(nodes/=nodes(j)))), &
                       save_pos=save_pos)
            else
              call set(upwind_values, i, nodes(j), maxval(field%val(nodes), mask=((nodes/=i).and.(nodes/=nodes(j)))), &
                       save_pos=save_pos)
            end if
            ! same for old field values
            if(node_val(old_field, i)<node_val(old_field, nodes(j))) then
              call set(old_upwind_values, i, nodes(j), minval(old_field%val(nodes), mask=((nodes/=i).and.(nodes/=nodes(j)))), &
                       save_pos=save_pos)
            else
              call set(old_upwind_values, i, nodes(j), maxval(old_field%val(nodes), mask=((nodes/=i).and.(nodes/=nodes(j)))), &
                       save_pos=save_pos)
            end if
         end do
      end do
      
   end subroutine calculate_upwind_values_local

  subroutine calculate_upwind_values_structured(x_field, field, upwind_values, &
                                           old_field, old_upwind_values)
      ! use the upwind value as close to the structured value as possible

      type(vector_field), intent(in) :: x_field  ! coordinates
      type(scalar_field), intent(in) :: field
      type(scalar_field), intent(in) :: old_field
      type(csr_matrix), intent(inout) :: upwind_values
      type(csr_matrix), intent(inout) :: old_upwind_values

      integer :: i, j, k
      integer, dimension(:), pointer :: nodes

      real, dimension(x_field%dim) :: vector1, vector2
      real :: dp, l_dp
      integer :: l_upwind_node

      integer :: save_pos=0 ! save the position in the matrix to optimise
      
      ewrite(2,*) 'in calculate_upwind_values_structured'
      ! the structured upwind value matrix is on the x_field mesh
      ! which cannot be periodic

      ! zero everything we have
      call zero(upwind_values)
      call zero(old_upwind_values)
      
      if(mesh_periodic(field)) then
        FLAbort("pseudo_structured_upwind_values don't work for periodic meshes")
      end if

      ! loop over nodes
      do i = 1, size(upwind_values, 1)
        ! find the neighbouring nodes using the matrix sparsity
        ! FIXME: the matrix sparsity is not the same as the mesh connectivity for 
        !        cube meshes!
        !        (this could potentially be very important for pseudo-structured upwind values
        !         as unnconected values may be used)
         nodes => row_m_ptr(upwind_values, i)

         ! loop over neighbouring nodes
         do j = 1, size(nodes)

            ! don't try using the donor as the downwind
            if(nodes(j)==i) cycle

            l_dp = 1.0  ! the maximum value it can take
            l_upwind_node = i  ! just in case we don't find a more suitable choice

            ! get the vector connecting the downwind and donor nodes
            vector1 = node_val(x_field, nodes(j))-node_val(x_field, i)
            vector1 = vector1/norm2(vector1) ! normalise

            do k = 1, size(nodes)

              ! don't try anything with the downwind or donor nodes themselves
              if((nodes(k)==nodes(j)).or.(nodes(k)==i)) cycle

              ! get the vector connecting the possible upwind and donor nodes
              vector2 = node_val(x_field, nodes(k))-node_val(x_field, i)
              vector2 = vector2/norm2(vector2) ! normalise

              ! take dot product
              dp = dot_product(vector1, vector2)

              ! we want to get the dot product that's as close to -1 as possible
              if (dp < l_dp) then
                l_dp = dp
                l_upwind_node = nodes(k)
              end if

            end do
            
            ! set the upwind value in the matrix
            call set(upwind_values, i, nodes(j), node_val(field, l_upwind_node), save_pos=save_pos)
            ! same for old field values
            call set(old_upwind_values, i, nodes(j), node_val(old_field, l_upwind_node), save_pos=save_pos)

         end do
      end do
      
   end subroutine calculate_upwind_values_structured

  subroutine calculate_boundary_normals(mesh, x, &
                                        normals, on_boundary, &
                                        surface_ids)

    ! calculate the outward pointing normals from the domain boundary
    ! and as a bonus tells you which nodes are on the boundary

    type(mesh_type), intent(inout) :: mesh
    type(vector_field), intent(in) :: x ! coordinates
    type(vector_field), intent(inout) :: normals
    logical, dimension(:), intent(inout) :: on_boundary
    ! an optional argument to specify on which surface ids the normal should be calculated
    integer, dimension(:), intent(in), optional :: surface_ids

    integer :: ele, sele, iloc, k
    integer, dimension(:), allocatable :: nodes_bdy
    real, dimension(:,:), allocatable :: x_ele_bdy, x_ele, normal_bdy

    type(cv_faces_type) :: cvfaces
    type(quadrature_type) :: bdyquad
    type(element_type) :: bdyshape

    logical, dimension(:), allocatable :: field_bc_type

    real :: normnormal

    if(mesh%shape%degree==0) then
      bdyquad = make_quadrature(vertices=face_vertices(mesh, 1), dim=(mesh_dim(mesh)-1), &
                                degree=1)
      bdyshape = make_element_shape(vertices=face_vertices(mesh,1), dim=(mesh_dim(mesh)-1), &
                                    degree=x%mesh%faces%shape%degree, quad=bdyquad)
    else
      cvfaces=find_cv_faces(vertices=ele_vertices(mesh,1), &
                            dimension=mesh_dim(mesh), &
                            polydegree=mesh%shape%degree, &
                            quaddegree=1)
      bdyshape=make_cvbdy_element_shape(cvfaces, x%mesh%faces%shape%degree)
    end if
    
    allocate(normal_bdy(mesh_dim(mesh), bdyshape%ngi))
    allocate(x_ele(mesh_dim(x),ele_loc(x,1)), &
             x_ele_bdy(mesh_dim(x),face_loc(x,1)))
    allocate(nodes_bdy(face_loc(mesh,1)))

    ! temporary hack
    assert(size(nodes_bdy)==size(normal_bdy, 2))  ! check that the no. of gauss pts. == no. of nodes
                                                  ! this should be guaranteed by quaddegree=1
                                                  ! but might break down for some elements

    allocate(field_bc_type(surface_element_count(mesh)))
    call get_periodic_boundary_condition(mesh, field_bc_type)

    ! calculate the normals at the nodes
    do sele = 1, surface_element_count(mesh)
    
      if(field_bc_type(sele)) cycle
      
      if(present(surface_ids)) then
        if(.not.any(surface_ids==surface_element_id(mesh, sele))) cycle
      end if

      ele = face_ele(x, sele)
      x_ele = ele_val(x, ele)
      x_ele_bdy = face_val(x, sele)

      nodes_bdy=face_global_nodes(normals, sele)

      on_boundary(nodes_bdy) = .true.

      if(mesh%shape%degree==0) then
        ! work out the normals at the gauss points of the face
        call transform_facet_to_physical(x, sele, normal=normal_bdy)
      else
        ! work out the normals at the gauss points of the face
        call transform_cvsurf_facet_to_physical(x_ele, x_ele_bdy, &
                              bdyshape, normal_bdy)
      end if

      ! here's where the temporary hack happens...
      ! add the gauss pt. normals to the nodes (they need to be the same length)
      call addto(normals, nodes_bdy, normal_bdy)

    end do ! sele

    ! normalise the normals
    do iloc = 1, node_count(normals)
    
      if(.not.on_boundary(iloc)) cycle

      normnormal = 0.0
      do k = 1, mesh_dim(normals)
        normnormal = normnormal + node_val(normals, k, iloc)**2
      end do
      normnormal = sqrt(normnormal)

      call set(normals, iloc, node_val(normals,iloc)/max(normnormal, tolerance))

    end do ! iloc

    if(mesh%shape%degree==0) then
      call deallocate(bdyquad)
    else
      call deallocate(cvfaces)
    end if
    call deallocate(bdyshape)

  end subroutine calculate_boundary_normals

  pure function project_upwind(i, j, x_field, &
                          on_boundary, normals) result(xc)

    ! for a given node pair calculate the vector connecting to the coordinates of a point
    ! just upwind

    integer, intent(in) :: i, j
    type(vector_field), intent(in) :: x_field
    logical, dimension(:), intent(in) :: on_boundary
    type(vector_field), intent(in) :: normals

    real, dimension(x_field%dim) :: xc

    integer :: k
    real, dimension(3,3) :: t_matrix, t_matrix_T
    real, dimension(3) :: n, vx, t2, t1, ref_vx

    if(on_boundary(i)) then  ! this is only true if you're on the boundary and have
                             ! requested that upwind values be reflected off it
       ! fill in zeros in case this is 2d
       n = 0.0
       vx = 0.0
       do k = 1, x_field%dim
         n(k) = normals%val(k)%ptr(i) ! extract the normal
         vx(k) = (x_field%val(k)%ptr(i) &
                 -x_field%val(k)%ptr(j))    ! extract the vector between the nodes
       end do
       t2 = cross_product(n, vx)      ! find the perpendicular surface tangent
       if(norm2(t2)<(1.e-5)*norm2(vx)) then   ! if the node pair is almost perpendicular to the surface
         do k = 1, x_field%dim                      ! then just reflect straight back
            xc(k) = c_distance*(x_field%val(k)%ptr(j) &
                                       -x_field%val(k)%ptr(i))
         end do
       else
         t2 = t2/norm2(t2)          ! normalise the surface tangent
         t1 = cross_product(n, -t2) ! get the orthogonal surface tangent
         ! you now have a coordinate system on the surface so
         ! create the reflected transformation matrix
         t_matrix(1,:) = -n ! the minus reflects the coordinate system back into the mesh
         t_matrix(2,:) = t1
         t_matrix(3,:) = t2
         ! create the transposed unreflected transformation matrix
         t_matrix_T(:,1) = n
         t_matrix_T(:,2) = t1
         t_matrix_T(:,3) = t2

         ! transform into the reflect coordinate system then immediately transform back to the original
         ! coordinates
         ref_vx=matmul(t_matrix_T, matmul(t_matrix, vx))
         do k = 1, x_field%dim
            ! find the vector to upwind coordinates
            xc(k) = c_distance*(ref_vx(k))
         end do
       end if
    else                    ! else just find the point outside the mesh (if you're actually on the boundary)
                            ! then (later) you'll just find the nearest element to this
                            ! and use the values from that
       do k = 1, x_field%dim
         ! find the vector to upwind coordinates
         xc(k) = -c_distance*(x_field%val(k)%ptr(j) &
                             -x_field%val(k)%ptr(i))
       end do
    end if

  end function project_upwind

  pure function project_vector(i, j, x, &
                          on_boundary, normals) result(d)

    ! for a given node pair calculate the direction of an upwind point

    integer, intent(in) :: i, j
    type(vector_field), intent(in) :: x
    logical, dimension(:), intent(in) :: on_boundary
    type(vector_field), intent(in) :: normals

    real, dimension(x%dim) :: d ! d is the vector between the donor and downwind nodes

    integer :: k
    real, dimension(3,3) :: t_matrix, t_matrix_T
    real, dimension(3) :: n, vx, t2, t1, ref_vx

    if(on_boundary(i)) then  ! this is only true if you're on the boundary and have
                             ! requested that upwind values be reflected off it
       ! fill in zeros in case this is 2d
       n = 0.0
       vx = 0.0
       do k = 1, x%dim
         n(k) = normals%val(k)%ptr(i) ! extract the normal
         vx(k) = (x%val(k)%ptr(i) &
                 -x%val(k)%ptr(j))    ! extract the vector between the nodes
       end do
       t2 = cross_product(n, vx)      ! find the perpendicular surface tangent
       if(norm2(t2)<(1.e-5)*norm2(vx)) then   ! if the node pair is almost perpendicular to the surface
         do k = 1, x%dim                      ! then just reflect straight back
            d(k) = x%val(k)%ptr(j)-x%val(k)%ptr(i)
         end do
       else
         t2 = t2/norm2(t2)          ! normalise the surface tangent
         t1 = cross_product(n, -t2) ! get the orthogonal surface tangent
         ! you now have a coordinate system on the surface so
         ! create the reflected transformation matrix
         t_matrix(1,:) = -n ! the minus reflects the coordinate system back into the mesh
         t_matrix(2,:) = t1
         t_matrix(3,:) = t2
         ! create the transposed unreflected transformation matrix
         t_matrix_T(:,1) = n
         t_matrix_T(:,2) = t1
         t_matrix_T(:,3) = t2

         ! transform into the reflect coordinate system then immediately transform back to the original
         ! coordinates
         ref_vx=matmul(t_matrix_T, matmul(t_matrix, vx))
         do k = 1, x%dim
            d(k) = ref_vx(k)
         end do
       end if
    else                    ! else just find the point outside the mesh (if you're actually on the boundary)
                            ! then (later) you'll just project to an upwind value outside... bounding won't be possible
      do k = 1, x%dim
         d(k) = x%val(k)%ptr(j)-x%val(k)%ptr(i)
      end do
    end if

  end function project_vector

  pure function calculate_volume_coordinates(x_ele, xc) result(coords)

    ! calculate the volume coordinates for simplex elements (tetrahedra)
    ! this is useful as they're the same as linear quadrature

    real, dimension(:,:), intent(in) :: x_ele
    real, dimension(:), intent(in) :: xc

    real, dimension(4) :: coords

    real, dimension(3) :: x1, x2, x3, x4
    real :: vol

    x1 = x_ele(:,1)
    x2 = x_ele(:,2)
    x3 = x_ele(:,3)
    x4 = x_ele(:,4)

    vol = calculate_simplex_volume()

    x1 = xc

    coords(1) = calculate_simplex_volume()

    x1 = x_ele(:,1)
    x2 = xc

    coords(2) = calculate_simplex_volume()

    x2 = x_ele(:,2)
    x3 = xc

    coords(3) = calculate_simplex_volume()

    x3 = x_ele(:,3)
    x4 = xc

    coords(4) = calculate_simplex_volume()

    coords = coords/vol

  contains

    pure function calculate_simplex_volume() result(svol)

      ! for a given set of nodes, calculate the volume

      real :: svol

      real, dimension(3) :: vector1, vector2, vector3

      vector1 = x2-x1
      vector2 = x3-x1
      vector3 = x4-x1

      svol = scalar_triple_product(vector1, vector2, vector3)
      svol = abs(svol)
!       svol = svol/6.0  ! this would make the volume correct
                         ! but as it gets cancelled out later
                         ! this is ommitted


    end function calculate_simplex_volume

  end function calculate_volume_coordinates

  function calculate_area_coordinates(x_ele, xc) result(coords)

    ! calculate the area coordinates for simplex elements (triangles)
    ! this is useful as they're the same as linear quadrature

    real, dimension(:,:), intent(in) :: x_ele
    real, dimension(:), intent(in) :: xc

    real, dimension(3) :: coords

    real, dimension(2) :: x1, x2, x3
    real :: area

    x1 = x_ele(:,1)
    x2 = x_ele(:,2)
    x3 = x_ele(:,3)

    area = calculate_simplex_area()

    x1 = xc

    coords(1) = calculate_simplex_area()

    x1 = x_ele(:,1)
    x2 = xc

    coords(2) = calculate_simplex_area()

    x2 = x_ele(:,2)
    x3 = xc

    coords(3) = calculate_simplex_area()

    coords = coords/area

    contains

    pure function calculate_simplex_area() result(sarea)

      ! for a given set of nodes, calculate the volume

      real :: sarea

      real, dimension(3) :: vector1, vector2, vector3

      vector1 = 0.0
      vector2 = 0.0

      vector1(1:2) = x2-x1
      vector2(1:2) = x3-x1

      vector3 = cross_product(vector1, vector2)
      sarea = norm2(vector3)
!       sarea = 0.5*sarea           ! this would make the area correct
                                  ! but as it gets cancelled out later
                                  ! this is ommitted

    end function calculate_simplex_area

  end function calculate_area_coordinates

  subroutine couple_upwind_values(upwind_values, old_upwind_values, cv_options)
  
    type(csr_matrix), intent(inout), dimension(:) :: upwind_values
    type(csr_matrix), intent(inout), dimension(:) :: old_upwind_values
    type(cv_options_type), dimension(:), intent(in) :: cv_options ! a wrapper type to pass in all the options for control volumes
    
    integer :: i, j, f, nfields, save_pos
    integer, dimension(:), pointer :: nodes
    real, dimension(size(upwind_values)) :: vals, old_vals
    
    nfields = size(upwind_values)
    save_pos = 0
    
    if(nfields>1) then
      do i = 1, size(upwind_values(1), 1)
        nodes => row_m_ptr(upwind_values(1), i)
        do j = 1, size(nodes)
          
          do f = 1, nfields
            vals(f) = val(upwind_values(f), i, nodes(j), save_pos=save_pos)
            old_vals(f) = val(old_upwind_values(f), i, nodes(j), save_pos=save_pos)
          end do
          
          do f = 2, nfields
            
            if (sum(vals(1:f))>cv_options(f)%sum_target_max) then
              vals(f) = cv_options(f)%sum_target_max-sum(vals(1:f-1))
              call set(upwind_values(f), i, nodes(j), vals(f), save_pos=save_pos)
            else if (sum(vals(1:f))<cv_options(f)%sum_target_min) then
              vals(f) = cv_options(f)%sum_target_min-sum(vals(1:f-1))
              call set(upwind_values(f), i, nodes(j), vals(f), save_pos=save_pos)
            end if
            
            if (sum(old_vals(1:f))>cv_options(f)%sum_target_max) then
              old_vals(f) = cv_options(f)%sum_target_max-sum(old_vals(1:f-1))
              call set(old_upwind_values(f), i, nodes(j), old_vals(f), save_pos=save_pos)
            else if (sum(old_vals(1:f))<cv_options(f)%sum_target_min) then
              old_vals(f) = cv_options(f)%sum_target_min-sum(old_vals(1:f-1))
              call set(old_upwind_values(f), i, nodes(j), old_vals(f), save_pos=save_pos)
            end if
          
          end do
        end do
      end do
    end if
  
  end subroutine couple_upwind_values

end module cv_upwind_values
