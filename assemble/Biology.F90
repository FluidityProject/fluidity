!    Copyright (C) 2008 Imperial College London and others.
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

module biology
  !!< This module implements a simple PZND (phytoplankton, zooplankton,
  !!< nutrient, detritus) model in ICOM.
  use fldebug
  use spud
  use global_parameters, only:FIELD_NAME_LEN,OPTION_PATH_LEN, PYTHON_FUNC_LEN
  use sparse_tools
  use fetools
  use fields
  use state_module
  use boundary_conditions
  use solvers
  use python_state
  use sparsity_patterns_meshes
  use field_options
  use fefields

  implicit none

  private
  public calculate_biology_terms, biology_check_options

  character(len=FIELD_NAME_LEN), dimension(6), parameter ::&
       & biology_fields = (/ &
       "Phytoplankton           ", &
       "Zooplankton             ", &
       "Nutrient                ", &
       "Detritus                ", &
       "Chlorophyll             ", &
       "Ammonium                "/)

  ! Boundary condition types:
  ! (the numbers should match up with the order in the 
  !  get_entire_boundary_condition call)
  integer :: BCTYPE_WEAKDIRICHLET=1

contains 

  subroutine calculate_biology_terms(state)
    !!< Calculate the incident light and hence the sources and sinks in the
    !!< biology equations.
    type(state_type), intent(inout) :: state

    character(len=OPTION_PATH_LEN) :: prefix, algorithm
    ! This is the photosynthetic radiation projected onto the 
    ! same mesh as the biology fields
    ! It also takes into account the "active" part of the solar radiation
    type(scalar_field) :: par_bio
    ! we use the phytoplankton as the "bio" mesh
    type(scalar_field), pointer :: phytoplankton, PhotosyntheticRadiation
    type(vector_field) :: coords
    integer :: stat

    call backup_source_terms(state)

    ! Don't do biology if it's not included in the model!
    if (.not.have_option("/ocean_biology")) return
    if (have_option("/ocean_biology/pznd")) then
       prefix="/ocean_biology/pznd"
       algorithm="/source_and_sink_algorithm"
    else if (have_option("/ocean_biology/lagrangian_ensemble")) then
       prefix="/ocean_biology/lagrangian_ensemble"
       algorithm="/biology_algorithm"
    else if (have_option("/ocean_biology/six_component")) then
       prefix="/ocean_biology/six_component"
       algorithm="/source_and_sink_algorithm"
    else
       FLExit("Unknown biology algorithm")
    end if

    ewrite(1,*) "Solving biology sources"

    ewrite(2,*) "will use ",trim(prefix)," model"

    ! Calculate the light field at every point.
    call solve_light_equation(state, prefix)

    par_bio = extract_scalar_field(state, "_PAR", stat)
    phytoplankton => extract_scalar_field(state, "Phytoplankton")
    if (stat /= 0) then
        ! field does not yet exist: create it
        call allocate(par_bio,phytoplankton%mesh, name="_PAR")
        call zero(par_bio)
        call insert(state, par_bio, par_bio%name)
        call deallocate(par_bio)
        par_bio = extract_scalar_field(state, "_PAR", stat)
    end if
    PhotosyntheticRadiation => extract_scalar_field(state, "PhotosyntheticRadiation")
    coords = get_coordinate_field(state, par_bio%mesh)
    ! project the Photosynthetic radaition field onto the _PAR field
    call project_field(PhotosyntheticRadiation, par_bio, coords)
    ! scale it to get the active part
    call scale(par_bio, 0.43)
    call deallocate(coords)

    ! Calculate the sources and sinks at every point.
    call calculate_biology_from_python(state, prefix, algorithm)

  end subroutine calculate_biology_terms


  subroutine calculate_biology_from_python(state, prefix, algorithm)
    ! Set the biological sources and sinks from python.
    type(state_type),intent(inout) :: state
    character(len=*), intent(in) :: prefix, algorithm
    character(len=PYTHON_FUNC_LEN) :: pycode

    if (.not.have_option(trim(prefix)//trim(algorithm))) then
       ! No sources and sinks specified.
       return
    end if
    
    ! Clean up to make sure that nothing else interferes
#ifdef HAVE_NUMPY
    call python_reset()
    call python_add_state(state)
    call get_option(trim(prefix)//trim(algorithm),pycode)
    ! And finally run the user's code
    call python_run_string(trim(pycode))    
#else
    ewrite(-1,*) "When configuring, make sure NumPy is found"
    FLExit("Python biological models require NumPy")
#endif

  end subroutine calculate_biology_from_python

  subroutine backup_source_terms(state)
    !!< Produce a backup copy of prescribed source fields and zero any
    !!< diagnostic source fields.
    type(state_type), intent(inout) :: state
    
    type(scalar_field) :: this_field, old_field
    integer :: i, stat
    
    field_loop: do i = 1, size(biology_fields)
       this_field= extract_scalar_field(state, &
            trim(biology_fields(i))//"Source", stat=stat)
       if (stat/=0) then
          cycle field_loop
       end if
       if (have_option(trim(this_field%option_path)//"/prescribed")) then
          ! Prescribed field.
          if (.not.has_scalar_field(state, &
               "Old"//trim(biology_fields(i))//"Source")) then

             call allocate(old_field, this_field%mesh, &
                  "Old"//trim(this_field%name))
             
             call set(old_field, this_field)
             call insert(state, old_field, trim(old_field%name))
             call deallocate(old_field)
          
          end if
       else
          ! Diagnostic field.
          
          call zero(this_field)
       end if
    end do field_loop

  end subroutine backup_source_terms

  subroutine solve_light_equation(state, prefix)
    !!< The photosynthetically active radiation at any depth in the ocean
    !!< is given by the equation: 
    !!<
    !!<  dL
    !!<  -- = (k_w - k_c P) L
    !!<  dg
    !!<
    !!< Where g is the direction of gravity, k is a constant and P is the
    !!< phytoplankton concentration.
    type(state_type), intent(inout) :: state
    !! Prefix to the options path. This selects the right biology model.
    character(len=*), intent(in) :: prefix
    !! Position
    type(vector_field) :: X
    !! Direction of gravity
    type(vector_field) :: g
    !! Light intensity
    type(scalar_field) :: light
    !! Phytoplankton density.
    type(scalar_field) :: P

    !! Normal first order sparsity pattern
    type(csr_sparsity), pointer :: sparsity
    !! Matrix for light
    type(csr_matrix) :: light_mat
    !! Right hand side for light equation
    type(scalar_field) :: rhs

    !! Field over the entire surface mesh containing bc values:
    type(scalar_field) :: bc_value
    !! Integer array of all surface elements indicating bc type
    !! (see below call to get_entire_boundary_condition):
    integer, dimension(:), allocatable :: bc_type
    
    ! Values of absorption coefficients for water and phytoplankton.
    real :: k_w, k_p
    
    ! Loop varable over elements
    integer :: ele

    ! If the PhotosyntheticRadiation field is not prognostic we don't solve
    ! anything.
    if (.not.have_option(trim(prefix)//&
         &"/scalar_field::PhotosyntheticRadiation/prognostic")) then
       return
    end if

    X=extract_vector_field(state, "Coordinate")
    g=extract_vector_field(state, "GravityDirection")
    light=extract_scalar_field(state, "PhotosyntheticRadiation")
    P=extract_scalar_field(state, "Phytoplankton")

    ! Only need first order sparsity as we have no diffusion
    sparsity=>get_csr_sparsity_firstorder(state, light%mesh, light%mesh)

    call get_option(trim(light%option_path)//&
         "/prognostic/absorption_coefficients/water", k_w)
    call get_option(trim(light%option_path)//&
         "/prognostic/absorption_coefficients/phytoplankton", k_p)

    call allocate(light_mat, sparsity)
    call zero(light_mat)

    call allocate(rhs, light%mesh, "LightRHS")
    call zero(rhs)
    
    ! Enquire about boundary conditions we're interested in
    ! Returns an integer array bc_type over the surface elements
    ! that indicates the bc type (in the order we specified, i.e.
    ! BCTYPE_WEAKDIRICHLET=1)
    allocate( bc_type(1:surface_element_count(light)) )
    call get_entire_boundary_condition(light, &
       & (/"weakdirichlet"/), &
       & bc_value, bc_type)

    assert(has_faces(X%mesh))
    assert(has_faces(light%mesh))

    do ele=1,element_count(light)
       call construct_light_element(ele, light_mat, rhs, X, g, P, light,&
            & bc_value, bc_type, k_w, k_p)
    end do

    call zero(light)
    call petsc_solve(light, light_mat, rhs)

    call deallocate(bc_value)
    call deallocate(light_mat)
    call deallocate(rhs)

  end subroutine solve_light_equation

  subroutine construct_light_element(ele, light_mat, rhs, X, g, P, light,&
       & bc_value, bc_type, k_w, k_p)
    integer, intent(in) :: ele
    type(csr_matrix), intent(inout) :: light_mat
    type(scalar_field), intent(inout) :: rhs
    type(vector_field), intent(in) :: X, g
    type(scalar_field), intent(in) :: P, light
    !! Field over the entire surface mesh containing bc values:
    type(scalar_field), intent(in):: bc_value
    !! Integer array of all surface elements indicating bc type
    !! (see above call to get_entire_boundary_condition):
    integer, dimension(:), intent(in):: bc_type    
    !! Values of absorption coefficients for water and phytoplankton.
    real :: k_w, k_p

    ! Neighbour element, face and neighbour face.
    integer :: ele_2, face, face_2
    ! Loops over faces.
    integer :: ni

    ! Variable transform times quadrature weights.
    real, dimension(ele_ngi(light,ele)) :: detwei
    ! Transformed gradient function for light.
    real, dimension(ele_loc(light, ele), ele_ngi(light, ele), &
         mesh_dim(light)) :: dl_t

    ! Gravity vector at each quadrature point.
    real, dimension(g%dim, ele_ngi(light, ele)) :: grav_q
    ! Phytoplankton concentration at each quadrature point.
    real, dimension(ele_ngi(P, ele)) :: P_q

    ! Neighbours of this element.
    integer, dimension(:), pointer :: neigh
    ! Whether the tracer field is continuous.
    logical :: dg

    ! Shape of the current element.
    type(element_type), pointer :: l_shape
    ! Global node numbers of the current element.
    integer, dimension(:), pointer :: light_ele

    dg=continuity(light)<0

    !----------------------------------------------------------------------
    ! Establish local node lists
    !----------------------------------------------------------------------

    light_ele=>ele_nodes(light,ele)  

    !----------------------------------------------------------------------
    ! Establish local shape functions
    !----------------------------------------------------------------------

    l_shape=>ele_shape(light, ele)

    ! Transform Tracer derivatives and weights into physical space.
    call transform_to_physical(X, ele, l_shape , dshape=dl_t, detwei=detwei)

    !----------------------------------------------------------------------
    ! Construct element-wise quantities.
    !----------------------------------------------------------------------

    grav_q = ele_val_at_quad(g, ele)

    P_q = k_w + k_p * max(ele_val_at_quad(P, ele), 0.0)

    !----------------------------------------------------------------------
    ! Local and global assembly in one fell swoop as it's so trivial.
    !----------------------------------------------------------------------

    call addto(light_mat, light_ele, light_ele, &
         ! dL
         ! --  (integrated by parts)
         ! dg
         -dshape_dot_vector_shape(dl_t, grav_q, l_shape, detwei)&
         ! RHS
         +shape_shape(l_shape, l_shape, detwei*P_q)&
         )

    !-------------------------------------------------------------------
    ! Interface integrals
    !-------------------------------------------------------------------
    
    neigh=>ele_neigh(light, ele)

    neighbourloop: do ni=1,size(neigh)

       !----------------------------------------------------------------------
       ! Find the relevant faces.
       !----------------------------------------------------------------------
       
       ! These finding routines are outside the inner loop so as to allow
       ! for local stack variables of the right size in
       ! construct_add_diff_interface_dg.

       ele_2=neigh(ni)
       
       ! Note that although face is calculated on field light, it is in fact
       ! applicable to any field which shares the same mesh topology.
       face=ele_face(light, ele, ele_2)
    
       if (ele_2>0) then
          ! Internal faces.

          if (.not. dg) then
             ! Continuous galerkin has no interior face terms.
             cycle neighbourloop
          else
             face_2=ele_face(light, ele_2, ele)
          end if
       else
          ! External face.
          face_2=face
       end if
       
       call construct_light_interface(face, face_2,&
            & light_mat, rhs, X, g, light, &
            & bc_value, bc_type)  
    end do neighbourloop

  end subroutine construct_light_element

  subroutine construct_light_interface(face, face_2,&
            & light_mat, rhs, X, g, light, bc_value, bc_type)
    !!< Construct the element boundary integrals on the ni-th face of
    !!< element ele. For continuous discretisation, this is only boundary
    !!< faces. For DG it's all of them.
    implicit none

    integer, intent(in) :: face, face_2
    type(csr_matrix), intent(inout) :: light_mat
    type(scalar_field), intent(inout) :: rhs
    ! We pass these additional fields to save on state lookups.
    type(vector_field), intent(in) :: X, g
    type(scalar_field), intent(in) :: light
    !! Field over the entire surface mesh containing bc values:
    type(scalar_field), intent(in):: bc_value
    !! Integer array of all surface elements indicating bc type
    !! (see above call to get_entire_boundary_condition):
    integer, dimension(:), intent(in):: bc_type

    ! Face objects and numberings.
    type(element_type), pointer ::l_shape, l_shape_2
    integer, dimension(face_loc(light,face)) :: l_face
    integer, dimension(face_loc(light,face_2)) :: l_face_2

    ! Note that both sides of the face can be assumed to have the same
    ! number of quadrature points.
    real, dimension(X%dim, face_ngi(X, face)) :: normal, grav_q
    real, dimension(face_ngi(X, face)) :: grav_flux
    logical, dimension(face_ngi(X, face)) :: influx
    ! Variable transform times quadrature weights.
    real, dimension(face_ngi(light,face)) :: detwei
    ! Gravity flux at each 
    ! Whether this is a boundary, and if so whether it is Dirichlet.
    logical :: boundary, Dirichlet

    ! Bilinear forms
    real, dimension(face_loc(light,face),face_loc(light,face)) ::&
         & nnlight_out 
    real, dimension(face_loc(light,face),face_loc(light,face_2)) ::&
         & nnlight_in 
    
    l_face=face_global_nodes(light, face)
    l_shape=>face_shape(light, face)

    l_face_2=face_global_nodes(light, face_2)
    l_shape_2=>face_shape(light, face_2)
    
    ! Boundary nodes have both faces the same.
    boundary=(face==face_2)
    dirichlet=.false.
    if (boundary .and. face < size(bc_type)) then
       if (bc_type(face)==BCTYPE_WEAKDIRICHLET) then
          dirichlet=.true.
       end if
    end if

    !----------------------------------------------------------------------
    ! Change of coordinates on face.
    !----------------------------------------------------------------------

    call transform_facet_to_physical(X, face, &
         &                          detwei_f=detwei,&
         &                          normal=normal) 
        
    !----------------------------------------------------------------------
    ! Construct element-wise quantities.
    !----------------------------------------------------------------------

    grav_q=face_val_at_quad(g, face)

    ! grav_flux is negative if gravity at this gauss point is directed
    ! into this element.
    grav_flux= sum(grav_q*normal,1)

    ! Flag for incoming flow.
    influx=grav_flux<0.0
   
    !----------------------------------------------------------------------
    ! Construct bilinear forms.
    !----------------------------------------------------------------------

    ! Calculate outflow boundary integral.
    nnlight_out=shape_shape(l_shape, l_shape,  &
         &                  merge(1.0,0.0,.not.influx)*grav_flux*detwei) 

    nnlight_in=shape_shape(l_shape, l_shape_2,  &
         &                  merge(1.0,0.0,influx)*grav_flux*detwei) 
    
    !----------------------------------------------------------------------
    ! Perform global assembly.
    !----------------------------------------------------------------------

    ! Insert flux terms in matrix.

    ! Outflux boundary integral.
    call addto(light_mat, l_face, l_face, nnlight_out)
       
    if (.not.dirichlet) then
       ! Influx boundary integral.
       call addto(light_mat, l_face, l_face_2,nnlight_in)
    end if

    ! Dirichlet boundary flux into rhs
       
    if (Dirichlet) then

       ! Inflow of Dirichlet value.
       call addto(RHS, l_face, &
            -matmul(nnlight_in,&
            ele_val(bc_value, face)))
 
    end if

  end subroutine construct_light_interface

  subroutine biology_check_options
    character(len=FIELD_NAME_LEN) :: buffer
    integer :: itmp, stat

    ! Don't do biology if it's not included in the model!
    if (.not.have_option("/ocean_biology/pznd") .or.  &
        .not.have_option("/ocean_biology/six_component")) return

    call get_option("/problem_type", buffer)
    if (buffer/="oceans") then
       FLExit("Biology modelling is only supported for problem type oceans.")
    end if

    if (.not.have_option("/physical_parameters/gravity")) then
       ewrite(-1, *) "Biology modelling requires gravity" 
       FLExit("(otherwise which way does the crap fall?)")
    end if

    if (have_option("/ocean_biology/pznd/scalar_field&
         &::PhotosyntheticRadiation/prognostic/solver/&
         &preconditioner::sor")) then
       ewrite(0, *) "Warning: Sor may not work for the PhotosyntheticRadiation "//&
            &"equation"
       ewrite(0, *) "Consider using ilu as a preconditioner instead."
    end if
    if (have_option("/ocean_biology/six_component/scalar_field&
         &::PhotosyntheticRadiation/prognostic/solver/&
         &preconditioner::sor")) then
       ewrite(0, *) "Warning: Sor may not work for the PhotosyntheticRadiation "//&
            &"equation"
       ewrite(0, *) "Consider using ilu as a preconditioner instead."
    end if

    call get_option("/timestepping/nonlinear_iterations", itmp, stat)
    
    if (stat/=0.or.itmp<2) then
       ewrite(0,*) "Warning: For stability reasons it is recommended that "//&
          "you have at least 2 nonlinear_iterations when using ocean biology"
    end if

  end subroutine biology_check_options

end module biology
