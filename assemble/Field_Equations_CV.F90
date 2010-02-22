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

module field_equations_cv
  !!< This module contains the assembly subroutines for advection
  !!< using control volumes
  use fields
  use sparse_matrices_fields
  use sparsity_patterns_meshes
  use state_module
  use spud
  use cv_shape_functions
  use cv_faces
  use cvtools
  use cv_fields
  use cv_upwind_values
  use cv_face_values
  use diagnostic_fields, only: calculate_diagnostic_variable
  use cv_options
  use diagnostic_variables, only: field_tag
  use boundary_conditions
  use boundary_conditions_from_options
  use divergence_matrix_cv, only: assemble_divergence_matrix_cv
  use global_parameters, only: OPTION_PATH_LEN
  use fefields, only: compute_lumped_mass
  use solvers, only: petsc_solve
  use transform_elements, only: transform_cvsurf_to_physical, &
                                transform_cvsurf_facet_to_physical
  use parallel_tools, only: getprocno
  use halos
  use field_options
  use state_fields_module

  implicit none

  private
  public :: solve_field_eqn_cv, field_equations_cv_check_options, &
            initialise_advection_convergence, calculate_auxiliary_gradient, coupled_cv_field_eqn, &
            assemble_advectiondiffusion_m_cv

  integer, dimension(:), allocatable, save :: conv_unit
  
  !! This allows a reference between the field and the file its meant to be
  !! writing to.
  character(len=OPTION_PATH_LEN), dimension(:), allocatable, save :: sfield_list

  ! are we moving the mesh?
  logical :: move_mesh

contains
    !************************************************************************
    ! solution wrapping subroutines
    subroutine solve_field_eqn_cv(field_name, state, global_it)
      !!< Construct and solve the advection-diffusion equation for the given
      !!< field using control volumes.

      !! Name of the field to be solved for.
      character(len=*), intent(in) :: field_name
      !! Collection of fields defining system state.
      type(state_type), dimension(:), intent(inout) :: state
      ! global iteration number - passed in so it can be output to file
      integer, intent(in) :: global_it

      ! Field to be solved for (plus its old and globally iterated versions)
      type(scalar_field), pointer :: tfield, oldtfield, it_tfield
      ! Density fields associated with tfield's equation (i.e. if its not pure advection)
      type(scalar_field), pointer :: tdensity, oldtdensity
      ! Coordinate field(s)
      type(vector_field), pointer :: x, x_old, x_new
      type(vector_field) :: x_tfield

      type(tensor_field), pointer :: diffusivity
      type(scalar_field), pointer :: source, absorption

      ! Change in tfield over one timestep.
      type(scalar_field) :: delta_tfield

      ! LHS equation matrix.
      type(csr_matrix) :: M
      ! Advection matrix
      type(csr_matrix) :: A_m
      ! Diffusion matrix
      type(csr_matrix) :: D_m
      ! sparsity structure to construct the matrices with
      type(csr_sparsity), pointer :: mesh_sparsity_1, mesh_sparsity, &
                                     mesh_sparsity_x, grad_m_t_sparsity

      ! Right hand side vector, lumped mass matrix, 
      ! locally iterated field (for advection iterations) 
      ! and local old field (for subcycling)
      type(scalar_field), pointer :: t_lumpedmass, q_lumpedmass
      type(scalar_field) :: t_lumpedmass_old, t_lumpedmass_new
      type(scalar_field) :: rhs, lumpedmass, advit_tfield, l_old_tfield
      ! Diffusion contribution to rhs
      type(scalar_field) :: diff_rhs

      ! local copy of option_path for solution field
      character(len=OPTION_PATH_LEN) :: option_path

      ! number of advection iterations and subcycles
      integer :: adv_iterations, no_subcycles
      ! iterators
      integer :: adv_it, sub
      ! time (to output to file), timestep, iterations tolerance, subcycling timestep
      real :: time, dt, error, adv_tolerance, sub_dt
      ! construct the matrix? construct the auxiliary diffusion eqn?
      logical :: getadvmat, getdiffmat, diffusion

      ! degree of quadrature to use on each control volume face
      integer :: quaddegree
      ! control volume face information
      type(cv_faces_type) :: cvfaces
      ! control volume shape function for volume and boundary
      type(element_type) :: u_cvshape, u_cvbdyshape
      type(element_type) :: ug_cvshape, ug_cvbdyshape
      type(element_type) :: x_cvshape, x_cvbdyshape, &
                            x_cvshape_full, x_cvbdyshape_full
      ! t_cvshape is the element with reduced numbers of derivatives
      ! taken across the control volume faces
      ! t_cvshape_full contains the derivatives with respect to the parent
      ! elements canonical coordinates evaluated at the control volume faces
      ! t_cvbdyshape_full is the same but on the boundary
      type(element_type) :: t_cvshape, t_cvshape_full, diff_cvshape_full, &
                                       t_cvbdyshape_full, diff_cvbdyshape_full

      ! options wrappers for tfield and tdensity
      type(cv_options_type) :: tfield_options, tdensity_options

      ! a dummy density in case we're solving for Advection
      type(scalar_field), pointer :: dummydensity, dummyscalar
      ! somewhere to put strings temporarily
      character(len=FIELD_NAME_LEN) :: tmpstring
      ! what equation type are we solving for?
      integer :: equation_type
      ! success indicators?
      integer :: stat, i
      ! the courant number field
      type(scalar_field) :: cfl_no
      ! nonlinear and grid velocities
      type(vector_field), pointer :: nu, ug
      ! advection velocity
      type(vector_field) :: advu
      !! Gravitational sinking term
      type(scalar_field), pointer :: sink
      !! Direction of gravity
      type(vector_field), pointer :: gravity

      ! assume explicitness?
      logical :: explicit
      ! if we're subcycling how fast can we go?
      real :: max_sub_cfl, max_cfl

      ! temporary hack to get around compiler failure to construct arrays of characters
      character(len=OPTION_PATH_LEN), dimension(1) :: option_path_array
      character(len=OPTION_PATH_LEN), dimension(1) :: density_option_path_array

      ewrite(2,*) 'in solve_field_eqn_cv'
      ewrite(2,*) 'solving for '//field_name//' in state '//trim(state(1)%name)

      ! extract lots of fields:
      ! the actual thing we're trying to solve for
      tfield=>extract_scalar_field(state(1), trim(field_name))
      ewrite_minmax(tfield%val)
      ! halo echange? - not currently necessary due to halo exchange directly after solve
      option_path=tfield%option_path
      ! its previous timelevel
      oldtfield=>extract_scalar_field(state(1), "Old"//trim(field_name))
      ! because fluidity resets tfield to oldtfield at the start of every
      ! global iteration we need to undo this so that the control volume faces
      ! are discretised using the most up to date values
      ! therefore extract the iterated values:
      it_tfield=>extract_scalar_field(state(1), "Iterated"//trim(field_name))
      ! and set tfield to them:
      call set(tfield, it_tfield)
      ewrite_minmax(tfield%val)
      ewrite_minmax(oldtfield%val)

      ! allocate dummy density in case density field isn't needed (this can be a constant field!)
      allocate(dummydensity)
      call allocate(dummydensity, tfield%mesh, name="DummyDensity", field_type=FIELD_TYPE_CONSTANT)
      call set(dummydensity, 1.0)
      dummydensity%option_path = " "

      ! find out equation type and hence if density is needed or not
      equation_type=equation_type_index(trim(option_path))
      select case(equation_type)
      case(FIELD_EQUATION_ADVECTIONDIFFUSION)
        ! density not needed so use a constant field for assembly
        tdensity=>dummydensity
        oldtdensity=>dummydensity
      case(FIELD_EQUATION_CONSERVATIONOFMASS, FIELD_EQUATION_REDUCEDCONSERVATIONOFMASS, &
           FIELD_EQUATION_INTERNALENERGY )
        call get_option(trim(option_path)//'/prognostic/equation[0]/density[0]/name', &
                        tmpstring)
        ! density needed so extract the type specified in the input
        ! ?? are there circumstances where this should be "Iterated"... need to be
        ! careful with priority ordering
        tdensity=>extract_scalar_field(state(1), trim(tmpstring))
        ewrite_minmax(tdensity%val)
        ! halo exchange? - not currently necessary when suboptimal halo exchange if density
        ! is solved for with this subroutine and the correct priority ordering.
        oldtdensity=>extract_scalar_field(state(1), "Old"//trim(tmpstring))
        ewrite_minmax(oldtdensity%val)
      end select

      ! now we can get the options for these fields
      ! handily wrapped in a new type...
      tfield_options=get_cv_options(tfield%option_path, tfield%mesh%shape%numbering%family)
      tdensity_options=get_cv_options(tdensity%option_path, tdensity%mesh%shape%numbering%family)

      ! extract fields from state
      nu=>extract_vector_field(state(1), "NonlinearVelocity")
      do i = 1, nu%dim
        ewrite_minmax(nu%val(i)%ptr)
      end do

      x=>extract_vector_field(state(1), "Coordinate")
      x_tfield=get_coordinate_field(state(1), tfield%mesh)

      ! find relative velocity
      call allocate(advu, nu%dim, nu%mesh, "AdvectionVelocity")
      call set(advu, nu)
      ! add in sinking velocity
      sink=>extract_scalar_field(state(1), trim(field_name)//"SinkingVelocity"&
           &, stat=stat)
      if(stat==0) then
        gravity=>extract_vector_field(state(1), "GravityDirection")
        ! this may perform a "remap" internally from CoordinateMesh to VelocitMesh
        call addto(advu, gravity, scale=sink)
      end if
      
      do i = 1, advu%dim
        ewrite_minmax(advu%val(i)%ptr)
      end do

      diffusion = .false.
      ! do we have a diffusivity - this will control whether we construct an auxiliary
      ! eqn or not (if BassiRebay is selected) or whether we construct gradients (if ElementGradient
      ! is selected)
      diffusivity=>extract_tensor_field(state(1), trim(field_name)//"Diffusivity", stat=stat)
      if(stat==0) diffusion = .true.

      ! allocate dummy scalar in case source/absorption fields aren't needed (this can be a constant field!)
      allocate(dummyscalar)
      call allocate(dummyscalar, tfield%mesh, name="DummyScalar", field_type=FIELD_TYPE_CONSTANT)
      call set(dummyscalar, 0.0)
      dummyscalar%option_path = " "

      source=>extract_scalar_field(state(1), trim(field_name)//"Source", stat=stat)
      if(stat/=0) then
        source=>dummyscalar
      else
        ewrite_minmax(source%val)
      end if
      absorption=>extract_scalar_field(state(1), trim(field_name)//"Absorption", stat=stat)
      if(stat/=0) then
        absorption=>dummyscalar
      else
        ewrite_minmax(absorption%val)
      end if

      ! create control volume shape functions
      call get_option("/geometry/quadrature/controlvolume_surface_degree", &
                     quaddegree, default=1)
      cvfaces=find_cv_faces(vertices=ele_vertices(tfield, 1), &
                            dimension=mesh_dim(tfield), &
                            polydegree=tfield%mesh%shape%degree, &
                            quaddegree=quaddegree)
      u_cvshape=make_cv_element_shape(cvfaces, nu%mesh%shape%degree)
      x_cvshape=make_cv_element_shape(cvfaces, x%mesh%shape%degree)
      t_cvshape=make_cv_element_shape(cvfaces, tfield%mesh%shape%degree)
      
      u_cvbdyshape=make_cvbdy_element_shape(cvfaces, nu%mesh%faces%shape%degree)
      x_cvbdyshape=make_cvbdy_element_shape(cvfaces, x%mesh%faces%shape%degree)

      if(diffusion.and.(tfield_options%diffusionscheme==CV_DIFFUSION_ELEMENTGRADIENT)) then
        x_cvshape_full=make_cv_element_shape(cvfaces, x%mesh%shape%degree, &
                                        type=ELEMENT_CONTROLVOLUME_SURFACE_BODYDERIVATIVES)
        t_cvshape_full=make_cv_element_shape(cvfaces, tfield%mesh%shape%degree, &
                                        type=ELEMENT_CONTROLVOLUME_SURFACE_BODYDERIVATIVES)
        diff_cvshape_full=make_cv_element_shape(cvfaces, diffusivity%mesh%shape%degree, &
                                        type=ELEMENT_CONTROLVOLUME_SURFACE_BODYDERIVATIVES)

        x_cvbdyshape_full=make_cvbdy_element_shape(cvfaces, x%mesh%faces%shape%degree, &
                                        type=ELEMENT_CONTROLVOLUME_SURFACE_BODYDERIVATIVES)
        t_cvbdyshape_full=make_cvbdy_element_shape(cvfaces, tfield%mesh%shape%degree, &
                                        type=ELEMENT_CONTROLVOLUME_SURFACE_BODYDERIVATIVES)
        diff_cvbdyshape_full=make_cvbdy_element_shape(cvfaces, diffusivity%mesh%shape%degree, &
                                        type=ELEMENT_CONTROLVOLUME_SURFACE_BODYDERIVATIVES)
      else
        x_cvshape_full=x_cvshape
        t_cvshape_full=t_cvshape
        diff_cvshape_full=t_cvshape
        x_cvbdyshape_full=x_cvbdyshape
        t_cvbdyshape_full=x_cvbdyshape
        diff_cvbdyshape_full=x_cvbdyshape

        call incref(x_cvshape_full)
        call incref(t_cvshape_full)
        call incref(diff_cvshape_full)
        call incref(x_cvbdyshape_full)
        call incref(t_cvbdyshape_full)
        call incref(diff_cvbdyshape_full)
      end if

      ! is this explicit?
      explicit=have_option(trim(option_path)//"/prognostic/explicit")
      ! find the timestep
      call get_option("/timestepping/timestep", dt)
      call get_option("/timestepping/current_time", time) ! so it can be output in the convergence file

      ! allocate and retrieve the cfl no. if necessary
      option_path_array(1) = trim(option_path)  ! temporary hack for
      density_option_path_array(1) = trim(tdensity%option_path) ! compiler failure
      call cv_disc_get_cfl_no(option_path_array, &
                      state(1), tfield%mesh, cfl_no, &
                      density_option_path_array)

      ! get the mesh sparsity for the matrices
      if(diffusion.and.(tfield_options%diffusionscheme==CV_DIFFUSION_BASSIREBAY)) then
        ! extend the sparsity
        mesh_sparsity => get_csr_sparsity_secondorder(state, tfield%mesh, diffusivity%mesh)

        mesh_sparsity_1 => get_csr_sparsity_firstorder(state, tfield%mesh, tfield%mesh)
        if(.not.(tfield%mesh==diffusivity%mesh)) then
          if(tfield%mesh%shape%degree>1) then
            FLAbort("To have a different diffusivity mesh the field must be at most P1")
          elseif(diffusivity%mesh%shape%degree>tfield%mesh%shape%degree) then
            FLAbort("The diffusivity mesh must be of a lower degree than the field")
          end if

          grad_m_t_sparsity => get_csr_sparsity_firstorder(state, tfield%mesh, diffusivity%mesh)
        else
          grad_m_t_sparsity => mesh_sparsity_1
        end if
        
      else

        mesh_sparsity => get_csr_sparsity_firstorder(state, tfield%mesh, tfield%mesh)

        grad_m_t_sparsity => mesh_sparsity

        mesh_sparsity_1 => mesh_sparsity
      end if

      if(mesh_periodic(tfield)) then
        if((tfield_options%upwind_scheme==CV_UPWINDVALUE_PROJECT_POINT).or.&
           (tfield_options%upwind_scheme==CV_UPWINDVALUE_PROJECT_GRAD)) then
           mesh_sparsity_x => get_csr_sparsity_firstorder(state, x_tfield%mesh, x_tfield%mesh)
        else
           mesh_sparsity_x => mesh_sparsity_1
        end if
      else
        mesh_sparsity_x => mesh_sparsity_1
      end if

      if(.not.explicit) then
        ! allocate the lhs matrix
        call allocate(M, mesh_sparsity, name=trim(field_name)//"Matrix")

        ! allocate the advection matrix
        call allocate(A_m, mesh_sparsity, name=trim(field_name)//"AdvectionMatrix")
        call zero(A_m)
      end if

      if(diffusion) then
        call allocate(D_m, sparsity=mesh_sparsity, name=trim(field_name)//"AuxiliaryMatrix")
        call zero(D_m)

        call allocate(diff_rhs, tfield%mesh, name=trim(field_name)//"DiffusionRHS")
        call zero(diff_rhs)
      end if

      ! allocate the rhs of the equation
      call allocate(rhs, tfield%mesh, name=trim(field_name)//"RHS")

      if(tfield%mesh%shape%degree>1) then
        ! try lumping on the submesh
        t_lumpedmass => get_lumped_mass_on_submesh(state, tfield%mesh)
      else
        ! then find the lumped mass
        ! (this replaces hart2/3d etc. in old code!!)
        t_lumpedmass => get_lumped_mass(state, tfield%mesh)
      end if
      ewrite_minmax(t_lumpedmass%val)

      move_mesh = have_option("/mesh_adaptivity/mesh_movement")
      if(move_mesh) then
        ewrite(2,*) "Moving mesh."
        x_old=>extract_vector_field(state(1), "OldCoordinate")
        x_new=>extract_vector_field(state(1), "IteratedCoordinate")
        call allocate(t_lumpedmass_old, tfield%mesh, name=trim(field_name)//"OldLumpedMass")
        call allocate(t_lumpedmass_new, tfield%mesh, name=trim(field_name)//"NewLumpedMass")
        if(tfield%mesh%shape%degree>1) then
          FLAbort("Lumping on submesh while moving the mesh not set up.")
        else
          call compute_lumped_mass(x_old, t_lumpedmass_old)
          call compute_lumped_mass(x_new, t_lumpedmass_new)
        end if
        ewrite_minmax(t_lumpedmass_old%val)
        ewrite_minmax(t_lumpedmass_new%val)
        
        ug=>extract_vector_field(state(1), "GridVelocity")
        do i = 1, ug%dim
          ewrite_minmax(ug%val(i)%ptr)
        end do

        ug_cvshape=make_cv_element_shape(cvfaces, ug%mesh%shape%degree)
        ug_cvbdyshape=make_cvbdy_element_shape(cvfaces, ug%mesh%faces%shape%degree)

      else
        ewrite(2,*) "Not moving mesh."
      end if

      call allocate(lumpedmass, tfield%mesh, name=trim(field_name)//"LumpedMass")
      call set(lumpedmass, t_lumpedmass)

      if(diffusion.and.(tfield_options%diffusionscheme==CV_DIFFUSION_BASSIREBAY)) then
        if(.not.(tfield%mesh==diffusivity%mesh)) then
          q_lumpedmass => get_lumped_mass(state, diffusivity%mesh)
        else
          q_lumpedmass => t_lumpedmass
        end if
      else
        q_lumpedmass => t_lumpedmass
      end if

      ! allocate a field to store the locally iterated values in
      call allocate(advit_tfield, tfield%mesh, name="AdvIterated"//trim(field_name))
      ! allocate a field to use as the local old field for subcycling
      call allocate(l_old_tfield, tfield%mesh, name="LocalOld"//trim(field_name))

      ! allocate a field to store the change between the old and new values
      call allocate(delta_tfield, tfield%mesh, name="Delta_"//trim(field_name))
      call zero(delta_tfield) ! Impose zero initial guess.
      ! Ensure delta_tfield inherits options from tfield for solver
      delta_tfield%option_path = option_path

      ! find out how many iterations we'll be doing
      call get_option(trim(option_path)//"/prognostic/temporal_discretisation&
                      &/control_volumes/number_advection_iterations", &
                      adv_iterations, default=1)
      call get_option(trim(option_path)//"/prognostic/temporal_discretisation&
                      &/control_volumes/number_advection_iterations/tolerance", &
                      adv_tolerance, default=0.0)

      sub_dt=dt  ! just in case I don't initialise this somehow
      ! are we subcycling?
      no_subcycles = 1
      call get_option(trim(option_path)//"/prognostic/temporal_discretisation&
                      &/control_volumes/number_advection_subcycles", &
                      no_subcycles, stat=stat)
      if(stat/=0) then
        ! have not specified a number of subcycles but perhaps we're using a 
        ! courant number definition?
        call get_option(trim(option_path)//"/prognostic/temporal_discretisation&
                        &/control_volumes/maximum_courant_number_per_subcycle", &
                        max_sub_cfl, stat=stat)
        if(stat==0) then
          max_cfl = maxval(cfl_no%val)
          call allmax(max_cfl)
          ! yes, we're subcycling
          ! we should have already calculated the courant number (or aborted in the attempt)
          no_subcycles=ceiling(max_cfl/max_sub_cfl)
          if(no_subcycles>1) then
            sub_dt=dt/real(no_subcycles)
            call scale(cfl_no, 1.0/real(no_subcycles))
          end if
        else
          ! no, we're not subcycling
          no_subcycles=1
          sub_dt = dt
        end if
      else
        if(no_subcycles>1) then
          sub_dt=dt/real(no_subcycles)
          call scale(cfl_no, 1.0/real(no_subcycles))
        end if
      end if

      ! when subcycling we're going to need to be starting each subcycle from the
      ! "new" old value but I don't want to screw with old code by updating the actual
      ! global timestep old value so lets create a copy now and update it instead
      call set(l_old_tfield, oldtfield)

      ewrite(2,*) 'entering subcycling_loop', no_subcycles
      ! subcycling loop
      subcycling_loop: do sub = 1, no_subcycles

        ! advection iteration loop
        advection_iteration_loop: do adv_it = 1, adv_iterations
          getadvmat=(adv_it==1).and.(sub==1).and.(.not.explicit)
          getdiffmat=(adv_it==1).and.(sub==1)

          ! record the value of tfield since the previous iteration
          call set(advit_tfield, tfield)

          if(move_mesh) then
            if(explicit) then
              call set(lumpedmass, t_lumpedmass_new)
            else
              call zero(M)
              call addto_diag(M, t_lumpedmass_new)
            end if
          else
            if(explicit) then
              call set(lumpedmass, t_lumpedmass)
            else
              call zero(M)
              call addto_diag(M, t_lumpedmass)
            end if
          end if
          call zero(rhs) ! this has to happen here rather than in the assembly routine
                         ! for back compatibility reasons

          ! assemble A_m and rhs
          call assemble_advectiondiffusion_m_cv(A_m, rhs, &
                                      tfield, l_old_tfield, tfield_options, &
                                      tdensity, oldtdensity, tdensity_options, &
                                      cvfaces, x_cvshape, x_cvbdyshape, &
                                      u_cvshape, u_cvbdyshape, t_cvshape, &
                                      ug_cvshape, ug_cvbdyshape, &
                                      state, advu, ug, x, x_tfield, cfl_no, &
                                      getadvmat, sub_dt, &
                                      mesh_sparsity_x, &
                                      diffusion=diffusion, getdiffmat=getdiffmat, &
                                      diffusivity=diffusivity, q_lumpedmass=q_lumpedmass, &
                                      D_m=D_m, diff_rhs=diff_rhs, grad_m_t_sparsity=grad_m_t_sparsity, &
                                      x_cvshape_full=x_cvshape_full, x_cvbdyshape_full=x_cvbdyshape_full, &
                                      t_cvshape_full=t_cvshape_full, diff_cvshape_full=diff_cvshape_full, &
                                      t_cvbdyshape_full=t_cvbdyshape_full, diff_cvbdyshape_full=diff_cvbdyshape_full)

          ! assemble it all into a coherent equation
          call assemble_field_eqn_cv(M, A_m, lumpedmass, rhs, &
                                    tfield, l_old_tfield, &
                                    tdensity, oldtdensity, &
                                    source, absorption, tfield_options%theta, &
                                    state, advu, sub_dt, explicit, &
                                    diffusion, D_m, diff_rhs, &
                                    t_lumpedmass_old, t_lumpedmass_new)


          ! Solve for the change in tfield.
          if(explicit) then
            call apply_dirichlet_conditions(lumpedmass, rhs, tfield, sub_dt)

            delta_tfield%val = rhs%val/lumpedmass%val
          else
            ! apply strong dirichlet boundary conditions (if any)
            ! note that weak conditions (known as control volume boundary conditions)
            ! will already have been applied
            call apply_dirichlet_conditions(M, rhs, tfield, sub_dt)

            call zero(delta_tfield)
            call petsc_solve(delta_tfield, M, rhs)
          end if

          ! reset tfield to l_old_tfield before applying change
          call set(tfield, l_old_tfield)
          ! Add the change in tfield to tfield.
          call addto(tfield, delta_tfield, sub_dt)

          call halo_update(tfield)  ! exchange the extended halos

          call test_and_write_advection_convergence(tfield, advit_tfield, x, &
                                    filename=trim(state(1)%name)//"__"//trim(tfield%name), &
                                    time=time+sub_dt, dt=sub_dt, it=global_it, adv_it=adv_it, &
                                    subcyc=sub, error=error)

          if(error<adv_tolerance) exit advection_iteration_loop

        end do advection_iteration_loop

        ! update the local old field to the new values and start again
        call set(l_old_tfield, tfield)

      end do subcycling_loop

      call deallocate(delta_tfield)
      call deallocate(advit_tfield)
      call deallocate(l_old_tfield)
      call deallocate(rhs)
      if(.not.explicit) call deallocate(A_m)
      if(.not.explicit) call deallocate(M)
      call deallocate(lumpedmass)
      call deallocate(cfl_no)
      call deallocate(x_cvbdyshape)
      call deallocate(x_cvbdyshape_full)
      call deallocate(u_cvbdyshape)
      call deallocate(x_cvshape)
      call deallocate(x_cvshape_full)
      call deallocate(u_cvshape)
      call deallocate(t_cvshape)
      call deallocate(t_cvshape_full)
      call deallocate(diff_cvshape_full)
      call deallocate(t_cvbdyshape_full)
      call deallocate(diff_cvbdyshape_full)
      call deallocate(cvfaces)
      call deallocate(advu)
      call deallocate(dummydensity)
      deallocate(dummydensity)
      call deallocate(dummyscalar)
      deallocate(dummyscalar)
      if (diffusion) then
        call deallocate(D_m)
        call deallocate(diff_rhs)
      end if
      call deallocate(x_tfield)
      if(move_mesh) then
        call deallocate(t_lumpedmass_new)
        call deallocate(t_lumpedmass_old)
        call deallocate(ug_cvshape)
        call deallocate(ug_cvbdyshape)
      end if

    end subroutine solve_field_eqn_cv

    subroutine coupled_cv_field_eqn(state, global_it)
      !!< This subroutine wraps the solve for groups of interdependent coupled fields.

      !! bucket full of fields from all materials
      type(state_type), dimension(:), intent(inout) :: state
      !! global iteration - passed in to be output to convergence file
      integer, intent(in) :: global_it

      type(scalar_field), pointer :: sfield

      !! list of fields that use the coupled_cv spatial_discretisation
      character(len=FIELD_NAME_LEN), dimension(:), allocatable :: field_name_list
      !! number of fields with the same name in different states that are therefore assumed to be interdependent
      integer, dimension(:), allocatable :: field_numbers

      integer :: i, f, c, nfields, nfield_groups

      ewrite(1,*) 'in coupled_cv_field_eqn'
      
      ! find the number of fields in all states (assumed to be in material_phases) that use the coupled_cv
      ! this is the maximum possible of fields we'll have to actually solve for (in reality will be fewer)
      nfields = option_count("/material_phase/scalar_field/prognostic/spatial_discretisation/coupled_cv")

      ewrite(2,*) 'nfields = ', nfields

      allocate(field_name_list(nfields))
      allocate(field_numbers(nfields))
      field_name_list = ""
      field_numbers = 0

      ! nfield_groups is the actual number of field groups we need to solve for
      nfield_groups = 0
      do i = 1, size(state)

        do f = 1, scalar_field_count(state(i))

          sfield => extract_scalar_field(state(i), f)
          if(have_option(trim(sfield%option_path)//"/prognostic/spatial_discretisation/coupled_cv")) then

            name_check_loop: do c = 1, nfields
              if(trim(sfield%name)==trim(field_name_list(c))) exit name_check_loop
            end do name_check_loop

            if(c>nfields) then
              ! not yet in the list so add it
              nfield_groups = nfield_groups + 1
              field_name_list(nfield_groups) = trim(sfield%name)
              field_numbers(nfield_groups) = 1
            else
              ! found in list, increment number of fields
              assert(c<=nfield_groups)
              field_numbers(c) = field_numbers(c) + 1
            end if

          end if

        end do

      end do
      
      ewrite(2,*) 'nfield_groups = ', nfield_groups

      do i = 1, nfield_groups
        assert(field_numbers(i)>0)
        call solve_coupled_cv(trim(field_name_list(i)), field_numbers(i), state, global_it)
      end do

      deallocate(field_name_list)
      deallocate(field_numbers)

    end subroutine coupled_cv_field_eqn

    subroutine solve_coupled_cv(field_name, nfields, state, global_it)
      !!< Construct and solve the advection equation for the given
      !!< field using coupled (i.e. interdependent face values) control volumes.

      !! Name of the field to be solved for.
      character(len=*), intent(in) :: field_name
      !! Number of interdependent fields (assumed to have the same name and be in different states)
      integer :: nfields
      !! Collection of fields defining system state.
      type(state_type), dimension(:), intent(inout) :: state
      !! global iteration number - passed in so it can be output to file
      integer, intent(in) :: global_it

      ! Field to be solved for (plus its old and globally iterated versions)
      type(scalar_field_pointer), dimension(nfields) :: tfield, oldtfield
      type(scalar_field), pointer :: it_tfield
      type(scalar_field), pointer :: tmpfield
      ! Density fields associated with tfield's equation (i.e. if its not pure advection)
      type(scalar_field_pointer), dimension(nfields) :: tdensity, oldtdensity
      ! Coordinate field
      type(vector_field), pointer :: x
      type(vector_field) :: x_tfield
      
      type(scalar_field_pointer), dimension(nfields) :: source, absorption

      ! Change in tfield over one timestep.
      type(scalar_field), dimension(nfields) :: delta_tfield

      ! LHS equation matrix.
      ! Advection matrix
      type(csr_matrix), dimension(nfields) :: M, A_m
      ! sparsity structure to construct the matrices with
      type(csr_sparsity), pointer :: mesh_sparsity, mesh_sparsity_x

      ! Right hand side vector, lumped mass matrix, 
      ! locally iterated field (for advection iterations) 
      ! and local old field (for subcycling)
      type(scalar_field), dimension(nfields) :: rhs, advit_tfield
      type(scalar_field_pointer), dimension(nfields) :: l_old_tfield
      type(scalar_field), dimension(nfields) :: lumpedmass
      type(scalar_field), pointer :: t_lumpedmass

      ! local copy of option_path for solution field
      character(len=OPTION_PATH_LEN), dimension(nfields) :: option_path
      character(len=OPTION_PATH_LEN), dimension(nfields) :: tdensity_option_path

      ! number of advection iterations and subcycles
      integer, dimension(nfields) :: adv_iterations, no_subcycles
      ! iterators
      integer :: adv_it, i, p, f, sub
      ! state indices
      integer, dimension(nfields) :: state_indices, tmp_state_indices, priorities
      ! time (to output to file), timestep, iterations tolerance, subcycling timestep
      real :: time, dt, sub_dt
      real, dimension(nfields) :: error
      real, dimension(nfields) :: adv_tolerance
      ! construct the matrix?
      logical, dimension(nfields) :: getmat

      ! degree of quadrature to use on each control volume face
      integer :: quaddegree
      ! control volume face information
      type(cv_faces_type) :: cvfaces
      ! control volume shape function for volume and boundary
      type(element_type) :: u_cvshape, u_cvbdyshape
      type(element_type) :: x_cvshape, x_cvbdyshape
      type(element_type) :: t_cvshape

      ! options wrappers for tfield and tdensity
      type(cv_options_type), dimension(nfields) :: tfield_options, tdensity_options

      ! a dummy density in case we're solving for Advection
      type(scalar_field), pointer :: dummydensity, dummyscalar
      ! somewhere to put strings temporarily
      character(len=FIELD_NAME_LEN) :: tmpstring
      ! what equation type are we solving for?
      integer :: equation_type
      ! success indicators?
      integer :: stat
      integer, dimension(nfields) :: cfl_sub_stat
      ! the courant number field
      type(scalar_field) :: cfl_no
      ! nonlinear and grid velocities
      type(vector_field), pointer :: nu
      ! advection velocity
      type(vector_field) :: advu
      ! assume explicitness?
      logical, dimension(nfields) :: explicit
      ! if we're subcycling how fast can we go?
      real, dimension(nfields) :: max_sub_cfl
      real :: max_cfl

      ewrite(1,*) 'in solve_coupled_cv'
      ewrite(2,*) 'solving for '//trim(field_name)//' '//int2str(nfields)//' times'

      ! find out where the fields are
      f = 0
      tmp_state_indices = 0
      priorities = 0
      do i = 1, size(state)
        tmpfield=>extract_scalar_field(state(i), trim(field_name))

        if(have_option(trim(tmpfield%option_path)//"/prognostic/spatial_discretisation/coupled_cv")) then
          f = f + 1
          call get_option(trim(tmpfield%option_path)//"/prognostic/priority", priorities(f))
          tmp_state_indices(f) = i
        end if
      end do

      assert(f==nfields)

      ! now work out the right order
      f = 0
      state_indices = 0
      do p = maxval(priorities), minval(priorities), -1
        do i=1, nfields
          if(priorities(i)==p) then
            f = f + 1
            state_indices(f) = tmp_state_indices(i)
          end if
        end do
      end do

      assert(f==nfields)

      ! allocate dummy density in case density field isn't needed (this can be a constant field!)
      allocate(dummydensity)
      call allocate(dummydensity, tmpfield%mesh, name="DummyDensity", field_type=FIELD_TYPE_CONSTANT)
      call set(dummydensity, 1.0)
      dummydensity%option_path = ""

      ! allocate dummy scalar in case source/absorption fields aren't needed (this can be a constant field!)
      allocate(dummyscalar)
      call allocate(dummyscalar, tmpfield%mesh, name="DummyScalar", field_type=FIELD_TYPE_CONSTANT)
      call zero(dummyscalar)
      dummyscalar%option_path = ""
      
      ! now extract everything in the right order
      do f = 1, nfields
        ewrite(2,*) 'extracting '//trim(field_name)//' from state '//trim(state(state_indices(f))%name)
        ! the field we want to solve for
        tfield(f)%ptr => extract_scalar_field(state(state_indices(f)), trim(field_name))
        ! its option path
        option_path(f)=tfield(f)%ptr%option_path
        ! its previous timelevel
        oldtfield(f)%ptr=>extract_scalar_field(state(state_indices(f)), "Old"//trim(field_name))
        ! because fluidity resets tfield to oldtfield at the start of every
        ! global iteration we need to undo this so that the control volume faces
        ! are discretised using the most up to date values
        ! therefore extract the iterated values:
        it_tfield=>extract_scalar_field(state(state_indices(f)), "Iterated"//trim(field_name))
        ! and set tfield to them:
        call set(tfield(f)%ptr, it_tfield)
        
        ! find out equation type and hence if density is needed or not
        equation_type=equation_type_index(trim(option_path(f)))
        select case(equation_type)
        case(FIELD_EQUATION_ADVECTIONDIFFUSION)
          ! density not needed so use a constant field for assembly
          tdensity(f)%ptr=>dummydensity
          oldtdensity(f)%ptr=>dummydensity
        case(FIELD_EQUATION_CONSERVATIONOFMASS, FIELD_EQUATION_REDUCEDCONSERVATIONOFMASS, &
            FIELD_EQUATION_INTERNALENERGY )
          call get_option(trim(option_path(f))//'/prognostic/equation[0]/density[0]/name', &
                          tmpstring)
          ! density needed so extract the type specified in the input
          ! ?? are there circumstances where this should be "Iterated"... need to be
          ! careful with priority ordering
          tdensity(f)%ptr=>extract_scalar_field(state(state_indices(f)), trim(tmpstring))
          ! halo exchange? - not currently necessary when suboptimal halo exchange if density
          ! is solved for with this subroutine and the correct priority ordering.
          oldtdensity(f)%ptr=>extract_scalar_field(state(state_indices(f)), "Old"//trim(tmpstring))
        end select
        ! its option path
        tdensity_option_path(f)=tdensity(f)%ptr%option_path

        ! now we can get the options for these fields
        ! handily wrapped in a new type...
        tfield_options(f)=get_cv_options(tfield(f)%ptr%option_path, tfield(f)%ptr%mesh%shape%numbering%family)
        tdensity_options(f)=get_cv_options(tdensity(f)%ptr%option_path, tdensity(f)%ptr%mesh%shape%numbering%family)

        source(f)%ptr=>extract_scalar_field(state(state_indices(f)), trim(field_name)//"Source", stat=stat)
        if(stat/=0) source(f)%ptr=>dummyscalar
        absorption(f)%ptr=>extract_scalar_field(state(state_indices(f)), trim(field_name)//"Absorption", stat=stat)
        if(stat/=0) absorption(f)%ptr=>dummyscalar

        ! is this explicit?
        explicit(f)=have_option(trim(option_path(f))//"/prognostic/explicit")

      end do
      
      ! we assume that all fields are on the same mesh as for the method to the work the faces must intersect!
      do f = 2, nfields
        assert(tfield(f)%ptr%mesh%shape%degree==tfield(1)%ptr%mesh%shape%degree)
        assert(all(tfield(f)%ptr%mesh%ndglno==tfield(1)%ptr%mesh%ndglno))
      end do
      
      ! for now we assume as this is a effectively a multimaterial problem 
      ! that all fields are advected using the same velocity
      ! on the same coordinate field
      ! extract velocity and coordinate fields from state
      nu=>extract_vector_field(state(state_indices(1)), "NonlinearVelocity")
      x=>extract_vector_field(state(state_indices(1)), "Coordinate")
      x_tfield = get_coordinate_field(state(state_indices(1)), tfield(1)%ptr%mesh)
      
      ! find relative velocity
      call allocate(advu, nu%dim, nu%mesh, "RelativeVelocity")
      call set(advu, nu)

      ! create control volume shape functions
      call get_option("/geometry/quadrature/controlvolume_surface_degree", &
                     quaddegree, default=1)
      cvfaces=find_cv_faces(vertices=ele_vertices(tfield(1)%ptr, 1), &
                            dimension=mesh_dim(tfield(1)%ptr), &
                            polydegree=element_degree(tfield(1)%ptr, 1), &
                            quaddegree=quaddegree)
      u_cvshape=make_cv_element_shape(cvfaces, nu%mesh%shape%degree)
      x_cvshape=make_cv_element_shape(cvfaces, x%mesh%shape%degree)
      t_cvshape=make_cv_element_shape(cvfaces, element_degree(tfield(1)%ptr, 1))
      u_cvbdyshape=make_cvbdy_element_shape(cvfaces, nu%mesh%faces%shape%degree)
      x_cvbdyshape=make_cvbdy_element_shape(cvfaces, x%mesh%faces%shape%degree)

      ! find the timestep
      call get_option("/timestepping/timestep", dt)
      call get_option("/timestepping/current_time", time) ! so it can be output in the convergence file

      ! allocate and retrieve the cfl no. if necessary
      call cv_disc_get_cfl_no(option_path, &
                      state(state_indices(1)), tfield(1)%ptr%mesh, cfl_no, &
                      tdensity_option_path)

      ! get the mesh sparsity for the matrices
      mesh_sparsity => get_csr_sparsity_firstorder(state, tfield(1)%ptr%mesh, tfield(1)%ptr%mesh)
      if(mesh_periodic(tfield(1)%ptr)) then
        if((tfield_options(1)%upwind_scheme==CV_UPWINDVALUE_PROJECT_POINT).or.&
           (tfield_options(1)%upwind_scheme==CV_UPWINDVALUE_PROJECT_GRAD)) then
          mesh_sparsity_x => get_csr_sparsity_firstorder(state, x_tfield%mesh, x_tfield%mesh)
        else
          mesh_sparsity_x => mesh_sparsity
        end if
      else
        mesh_sparsity_x => mesh_sparsity
      end if


      do f = 1, nfields
        ! allocate the lhs matrix
        if(.not.explicit(f)) then
          call allocate(M(f), mesh_sparsity, name=trim(field_name)//"Matrix")
          call zero(M(f))
  
          ! allocate the advection matrix
          call allocate(A_m(f), mesh_sparsity, name=trim(field_name)//int2str(f)//"AdvectionMatrix")
          call zero(A_m(f))
        end if

        ! allocate the rhs of the equation
        call allocate(rhs(f), tfield(f)%ptr%mesh, name=trim(field_name)//int2str(f)//"RHS")
      end do

      if(element_degree(tfield(1)%ptr, 1)>1) then
        ! try lumping on the submesh
        t_lumpedmass => get_lumped_mass_on_submesh(state, tfield(1)%ptr%mesh)
      else
        ! then find the lumped mass
        ! (this replaces hart2/3d etc. in old code!!)
        t_lumpedmass => get_lumped_mass(state, tfield(1)%ptr%mesh)
      end if
      ewrite_minmax(t_lumpedmass%val)

      move_mesh = have_option("/mesh_adaptivity/mesh_movement")
      if(move_mesh) then
        FLAbort("Moving meshes not set up with coupled cv.")
      end if

      do f = 1, nfields
        call allocate(lumpedmass(f), tfield(1)%ptr%mesh, name=trim(field_name)//"LocalLumpedMass")
        call set(lumpedmass(f), t_lumpedmass) ! save it in case we're explicit
      end do

      do f = 1, nfields
        ! allocate a field to store the locally iterated values in
        call allocate(advit_tfield(f), tfield(f)%ptr%mesh, name="AdvIterated"//int2str(f)//trim(field_name))
        ! allocate a field to use as the local old field for subcycling
        allocate(l_old_tfield(f)%ptr)
        call allocate(l_old_tfield(f)%ptr, tfield(f)%ptr%mesh, name="LocalOld"//int2str(f)//trim(field_name))
        ! when subcycling we're going to need to be starting each subcycle from the
        ! "new" old value but I don't want to screw with old code by updating the actual
        ! global timestep old value so lets create a copy now and update it instead
        call set(l_old_tfield(f)%ptr, oldtfield(f)%ptr)

        ! allocate a field to store the change between the old and new values
        call allocate(delta_tfield(f), tfield(f)%ptr%mesh, name="Delta_"//int2str(f)//trim(field_name))
        call zero(delta_tfield(f)) ! Impose zero initial guess.
        ! Ensure delta_tfield inherits options from tfield for solver
        delta_tfield(f)%option_path = option_path(f)
      end do

      adv_iterations = 1
      adv_tolerance = 0.0
      no_subcycles = 1
      cfl_sub_stat = 1
      sub_dt=dt  ! just in case I don't initialise this somehow
      do f = 1, nfields
        ! find out how many iterations we'll be doing
        call get_option(trim(option_path(f))//"/prognostic/temporal_discretisation&
                        &/control_volumes/number_advection_iterations", &
                        adv_iterations(f), default=1)

        call get_option(trim(option_path(f))//"/prognostic/temporal_discretisation&
                        &/control_volumes/number_advection_iterations/tolerance", &
                        adv_tolerance(f), default=0.0)

        call get_option(trim(option_path(f))//"/prognostic/temporal_discretisation&
                        &/control_volumes/number_advection_subcycles", &
                        no_subcycles(f), stat=cfl_sub_stat(f))

      end do
      assert(all(adv_iterations==adv_iterations(1)))
      assert(all(adv_tolerance==adv_tolerance(1)))
      assert(all(no_subcycles==no_subcycles(1)))
      assert(all(cfl_sub_stat==cfl_sub_stat(1)))

      stat = cfl_sub_stat(1)
      cfl_sub_stat = 1
      max_sub_cfl=0.0
      if(stat/=0) then
        ! have not specified a number of subcycles but perhaps we're using a 
        ! courant number definition?
        do f = 1, nfields
          call get_option(trim(option_path(f))//"/prognostic/temporal_discretisation&
                          &/control_volumes/maximum_courant_number_per_subcycle", &
                          max_sub_cfl(f), stat=cfl_sub_stat(f))
        end do
        assert(all(max_sub_cfl==max_sub_cfl(1)))
        assert(all(cfl_sub_stat==cfl_sub_stat(1)))
        if(cfl_sub_stat(1)==0) then
          max_cfl = maxval(cfl_no%val)
          call allmax(max_cfl)
          ! yes, we're subcycling
          ! we should have already calculated the courant number (or aborted in the attempt)
          no_subcycles=ceiling(max_cfl/max_sub_cfl(1))
          if(no_subcycles(1)>1) then
            sub_dt=dt/real(no_subcycles(1))
            call scale(cfl_no, 1.0/real(no_subcycles(1)))
          end if
        else
          ! no, we're not subcycling
          no_subcycles=1
          sub_dt = dt
        end if
      else
        if(no_subcycles(1)>1) then
          sub_dt=dt/real(no_subcycles(1))
          call scale(cfl_no, 1.0/real(no_subcycles(1)))
        end if
      end if

      ewrite(2,*) 'entering subcycling_loop', no_subcycles(1)
      ! subcycling loop
      subcycling_loop: do sub = 1, no_subcycles(1)

        ! advection iteration loop
        advection_iteration_loop: do adv_it = 1, adv_iterations(1)

          do f = 1, nfields
            getmat(f)=(adv_it==1).and.(sub==1).and.(.not.explicit(f))
            
            ! record the value of tfield since the previous iteration
            call set(advit_tfield(f), tfield(f)%ptr)

            if(explicit(f)) then
              call set(lumpedmass(f), t_lumpedmass)
            else
              call zero(M(f))
              call addto_diag(M(f), t_lumpedmass)
            end if
            call zero(rhs(f))
          end do

          ! assemble A_m and rhs
          call assemble_coupled_advection_m_cv(A_m, rhs, &
                                      tfield, l_old_tfield, tfield_options, &
                                      tdensity, oldtdensity, tdensity_options, &
                                      cvfaces, x_cvshape, x_cvbdyshape, &
                                      u_cvshape, u_cvbdyshape, t_cvshape, &
                                      state, advu, x, x_tfield, cfl_no, &
                                      getmat, sub_dt, &
                                      mesh_sparsity_x)



          do f = 1, nfields

            ! assemble it all into a coherent equation
            call assemble_field_eqn_cv(M(f), A_m(f), lumpedmass(f), rhs(f), &
                                      tfield(f)%ptr, l_old_tfield(f)%ptr, &
                                      tdensity(f)%ptr, oldtdensity(f)%ptr, &
                                      source(f)%ptr, absorption(f)%ptr, tfield_options(f)%theta, &
                                      state(state_indices(f):state_indices(f)), advu, sub_dt, explicit(f))

            ! Solve for the change in tfield.
            if(explicit(f)) then
              call apply_dirichlet_conditions(lumpedmass(f), rhs(f), tfield(f)%ptr, sub_dt)

              delta_tfield(f)%val = rhs(f)%val/lumpedmass(f)%val
            else
              ! apply strong dirichlet boundary conditions (if any)
              ! note that weak conditions (known as control volume boundary conditions)
              ! will already have been applied
              call apply_dirichlet_conditions(M(f), rhs(f), tfield(f)%ptr, sub_dt)

              call zero(delta_tfield(f))
              call petsc_solve(delta_tfield(f), M(f), rhs(f))
            end if

            ewrite_minmax(delta_tfield(f)%val)

            ! reset tfield to l_old_tfield before applying change
            call set(tfield(f)%ptr, l_old_tfield(f)%ptr)
            ! Add the change in tfield to tfield.
            call addto(tfield(f)%ptr, delta_tfield(f), sub_dt)

            call halo_update(tfield(f)%ptr)  ! exchange the extended halos

            call test_and_write_advection_convergence(tfield(f)%ptr, advit_tfield(f), x, &
                                      filename=trim(state(state_indices(f))%name)//"__"//trim(tfield(f)%ptr%name), &
                                      time=time+sub_dt, dt=sub_dt, it=global_it, adv_it=adv_it, &
                                      subcyc=sub, error=error(f))

          end do

          if(all(error<adv_tolerance)) exit advection_iteration_loop

        end do advection_iteration_loop

        do f = 1, nfields
          ewrite_minmax(tfield(f)%ptr%val)
          ! update the local old field to the new values and start again
          call set(l_old_tfield(f)%ptr, tfield(f)%ptr)
        end do

      end do subcycling_loop

      do f = 1, nfields
        call deallocate(delta_tfield(f))
        call deallocate(advit_tfield(f))
        call deallocate(l_old_tfield(f)%ptr)
        deallocate(l_old_tfield(f)%ptr)
        call deallocate(rhs(f))
        if(.not.explicit(f)) call deallocate(A_m(f))
        call deallocate(lumpedmass(f))
        if(.not.explicit(f)) call deallocate(M(f))
      end do
      call deallocate(cfl_no)
      call deallocate(x_cvbdyshape)
      call deallocate(u_cvbdyshape)
      call deallocate(x_cvshape)
      call deallocate(u_cvshape)
      call deallocate(t_cvshape)
      call deallocate(cvfaces)
      call deallocate(advu)
      call deallocate(dummydensity)
      deallocate(dummydensity)
      call deallocate(dummyscalar)
      deallocate(dummyscalar)
      call deallocate(x_tfield)

    end subroutine solve_coupled_cv
    ! end of solution wrapping subroutines
    !************************************************************************
    !************************************************************************
    ! equation wrapping subroutines
    subroutine assemble_field_eqn_cv(M, A_m, lumpedmass, rhs, &
                                    tfield, oldtfield, &
                                    tdensity, oldtdensity, &
                                    source, absorption, theta, &
                                    state, advu, dt, explicit, &
                                    diffusion, D_m, diff_rhs, &
                                    lumpedmass_old, lumpedmass_new)

      ! This subroutine assembles the equation
      ! M(T^{n+1}-T^{n})/\Delta T = rhs
      ! for control volumes.
      ! By the time you get here M should already contain the mass
      ! components (and if back compatible the diffusional components)
      ! of the equation.

      ! inputs/outputs:
      ! lhs matrix
      type(csr_matrix), intent(inout) :: M
      ! matrix containing advective terms - to be incorporated
      ! into M during this subroutine
      type(csr_matrix), intent(inout) :: A_m
      ! rhs of equation
      type(scalar_field), intent(inout) :: lumpedmass, rhs
      ! the field we are solving for
      type(scalar_field), intent(inout) :: tfield
      type(scalar_field), intent(inout) :: oldtfield, tdensity, oldtdensity
      type(scalar_field), intent(inout) :: source, absorption
      ! time discretisation parameter
      real, intent(in) :: theta
      ! bucket full of fields
      type(state_type), dimension(:), intent(inout) :: state
      ! advection velocity
      type(vector_field), intent(inout) :: advu
      ! the timestep
      real, intent(in) :: dt
      ! are we assuming this is a fully explicit equation?
      logical, intent(in) :: explicit
      ! diffusion:
      logical, intent(in), optional :: diffusion
      type(csr_matrix), intent(inout), optional :: D_m
      type(scalar_field), intent(inout), optional :: diff_rhs
      type(scalar_field), intent(in), optional :: lumpedmass_old
      type(scalar_field), intent(in), optional :: lumpedmass_new

      ! local memory:
      ! for all equation types:
      ! product of A_m or D_m and oldtfield
      type(scalar_field) :: MT_old
      type(scalar_field) :: masssource, massabsorption, massconservation

      ! for InternalEnergy equations:
      ! sparsity for CT_m
      type(csr_sparsity), pointer :: gradient_sparsity
      ! divergence matrix for energy equation
      type(block_csr_matrix) :: CT_m
      ! pressure
      type(scalar_field), pointer :: p
      ! the assembled pressure term
      type(scalar_field) :: pterm
      ! atmospheric pressure for the energy equation
      real :: atmospheric_pressure

      ! for ConservationOfMass equations:
      ! tmpfield used in assembly of conservation term, consterm
      type(scalar_field) :: consterm, tmpfield

      ! for ReducedConservationOfMass equations:
      real :: tdensity_theta

      ! self explanatory strings
      character(len=FIELD_NAME_LEN) :: tmpstring
      integer :: equation_type

      logical :: l_diffusion

      if(present(diffusion)) then
        l_diffusion = diffusion
      else
        l_diffusion = .false.
      end if

      if(l_diffusion) then
        if(.not.present(D_m).or..not.present(diff_rhs)) then
          FLAbort("Must supply a diffusion matrix and rhs to use diffusion.")
        end if
      end if

      ! allocate some memory for assembly
      call allocate(MT_old, rhs%mesh, name="MT_oldProduct" )
      call allocate(masssource, rhs%mesh, name="MassSourceProduct" )
      call set(masssource, lumpedmass)
      call scale(masssource, source)
      call allocate(massabsorption, rhs%mesh, name="MassAbsorptionProduct" )
      call set(massabsorption, lumpedmass)
      call scale(massabsorption, absorption)
      
      if(move_mesh) then
        if((.not.present(lumpedmass_new)).or.(.not.present(lumpedmass_old))) then
          FLAbort("I need a new and an old lumped mass to move the mesh.")
        end if
        call allocate(massconservation, rhs%mesh, name="MovingMeshMassConservation")
        call set(massconservation, lumpedmass_old)
        call addto(massconservation, lumpedmass_new, scale=-1.0)
        call scale(massconservation, 1./dt)
        call scale(massconservation, oldtfield)
      end if

      ! find out equation type and hence if density is needed or not
      equation_type=equation_type_index(trim(tfield%option_path))
      ! now we need to incorporate A_m into M and turn the equation into
      ! rate of change form (as well as adding in any extra terms for InternalEnergy
      ! for instance)
      select case(equation_type)
      case (FIELD_EQUATION_ADVECTIONDIFFUSION)

        ! [M + A_m](T^{n+1}-T^{n})/dt = rhs - A_m*T^{n}

        if(.not.explicit) then
          ! construct M
          call addto(M, A_m, dt)
          call addto_diag(M, massabsorption, theta*dt)
          
          ! construct rhs
          call mult(MT_old, A_m, oldtfield)
          call addto(rhs, MT_old, -1.0)
        end if


        call addto(rhs, masssource)

        ! massabsorption has already been added to the matrix so it can now be scaled
        ! by the old field value to add it to the rhs
        call scale(massabsorption, oldtfield)
        call addto(rhs, massabsorption, -1.0)

        if(l_diffusion) then
          call mult(MT_old, D_m, oldtfield)
          call addto(rhs, MT_old, -1.0)
          call addto(rhs, diff_rhs, -1.0)

          if(.not.explicit) then
            call addto(M, D_m, theta*dt)
          end if
        end if
        
        if(move_mesh) then
          call addto(rhs, massconservation)
        end if

      case (FIELD_EQUATION_CONSERVATIONOFMASS)

        ! [\rho^{n+1}M + dt*A_m](T^{n+1}-T^{n})/dt = rhs - A_m*T^{n} - M*(\rho^{n+1}-\rho^{n})*T^{n}/dt

        ! construct rhs
        call allocate(tmpfield, tfield%mesh, name="DensityDifference")
        call allocate(consterm, tdensity%mesh, name="DensityDifference")
        tmpfield%val=(1./dt)*(tdensity%val-oldtdensity%val)*oldtfield%val

        call allocate(consterm, lumpedmass%mesh, "ConservationTerm")
        call set(consterm, lumpedmass)
        call scale(consterm, tmpfield)
        call addto(rhs, consterm, -1.0)

        call deallocate(tmpfield)
        call deallocate(consterm)
        
        ! construct M:
        ! multiply the diagonal of M by the up to date density
        if(explicit) then
          lumpedmass%val = lumpedmass%val*tdensity%val
        else
          call mult_diag(M, tdensity)
          
          call addto(M, A_m, dt)
          call addto_diag(M, massabsorption, theta*dt)
        
          call mult(MT_old, A_m, oldtfield)
          call addto(rhs, MT_old, -1.0)
        end if
        
        call addto(rhs, masssource)

        ! massabsorption has already been added to the matrix so it can now be scaled
        ! by the old field value to add it to the rhs
        call scale(massabsorption, oldtfield)
        call addto(rhs, massabsorption, -1.0)
        
        if(move_mesh) then
          FLAbort("Moving mesh with this equation type not yet supported.")
        end if

      case (FIELD_EQUATION_REDUCEDCONSERVATIONOFMASS)

        ! [\rho^{n}M + dt*A_m](T^{n+1}-T^{n})/dt = rhs - A_m*T^{n}

        call get_option(trim(tdensity%option_path)//"/prognostic/temporal_discretisation&
                             &/theta", tdensity_theta)

        ! construct M
        ! multiply the diagonal by the previous timesteps density
        if(explicit) then
          lumpedmass%val = lumpedmass%val*(tdensity_theta*tdensity%val+(1.0-tdensity_theta)*oldtdensity%val)
        else
          call mult_diag(M, ((tdensity_theta)*tdensity%val+(1.0-tdensity_theta)*oldtdensity%val))
          call addto(M, A_m, dt)
          call addto_diag(M, massabsorption, theta*dt)
        
          ! construct rhs
          call mult(MT_old, A_m, oldtfield)
          call addto(rhs, MT_old, -1.0)
        end if

        call addto(rhs, masssource)

        ! massabsorption has already been added to the matrix so it can now be scaled
        ! by the old field value to add it to the rhs
        call scale(massabsorption, oldtfield)
        call addto(rhs, massabsorption, -1.0)

        if(move_mesh) then
          FLAbort("Moving mesh with this equation type not yet supported.")
        end if

      case (FIELD_EQUATION_INTERNALENERGY)

        ! [\rho^{n+1}M + dt*A_m](T^{n+1}-T^{n})/dt = rhs - A_m*T^{n} - (p+atm_p)*CT_m*u

        ! construct rhs
        p=>extract_scalar_field(state(1), "Pressure")
        ewrite_minmax(p%val)
        assert(p%mesh==tfield%mesh)
        ! halo exchange not necessary as it is done straight after solve
        call get_option(trim(p%option_path)//'/prognostic/atmospheric_pressure', &
                              atmospheric_pressure, default=0.0)
        gradient_sparsity => get_csr_sparsity_firstorder(state, p%mesh, advu%mesh)

        call allocate(CT_m, gradient_sparsity, (/1, advu%dim/), name="DivergenceMatrix" )
        call assemble_divergence_matrix_cv(CT_m, state(1), &
                                           test_mesh=p%mesh, field=advu)

        call allocate(pterm, p%mesh, "PressureTerm")

        ! construct the pressure term
        call mult(pterm, CT_m, advu) 
                                ! should this really be the advection velocity or just the relative or the nonlinear?
        pterm%val = pterm%val*(p%val+atmospheric_pressure)

        call addto(rhs, pterm, -1.0)

        call deallocate(CT_m)
        call deallocate(pterm)

        ! construct M
        if(explicit) then
          lumpedmass%val = lumpedmass%val*tdensity%val
        else
          call mult_diag(M, tdensity)
          call addto(M, A_m, dt)
          call addto_diag(M, massabsorption, theta*dt)
        
          call mult(MT_old, A_m, oldtfield)
          call addto(rhs, MT_old, -1.0)
        end if

        call addto(rhs, masssource)

        ! massabsorption has already been added to the matrix so it can now be scaled
        ! by the old field value to add it to the rhs
        call scale(massabsorption, oldtfield)
        call addto(rhs, massabsorption, -1.0)

        if(move_mesh) then
          FLAbort("Moving mesh with this equation type not yet supported.")
        end if

      end select

      call deallocate(masssource)
      call deallocate(massabsorption)
      call deallocate(MT_old)
      if(move_mesh) then
        call deallocate(massconservation)
      end if

    end subroutine assemble_field_eqn_cv
    ! end of equation wrapping subroutines
    !************************************************************************
    !************************************************************************
    ! assembly subroutines 
    subroutine assemble_advectiondiffusion_m_cv(A_m, rhs, &
                                       tfield, oldtfield, tfield_options, &
                                       tdensity, oldtdensity, tdensity_options, &
                                       cvfaces, x_cvshape, x_cvbdyshape, &
                                       u_cvshape, u_cvbdyshape, t_cvshape, &
                                       ug_cvshape, ug_cvbdyshape, &
                                       state, advu, ug, x, x_tfield, cfl_no, &
                                       getadvmat, dt, &
                                       mesh_sparsity, &
                                       diffusion, getdiffmat, diffusivity, q_lumpedmass, &
                                       D_m, diff_rhs, grad_m_t_sparsity, &
                                       x_cvshape_full, x_cvbdyshape_full, &
                                       t_cvshape_full, diff_cvshape_full, &
                                       t_cvbdyshape_full, diff_cvbdyshape_full, &
                                       reference_field)

      ! This subroutine assembles the advection matrix and rhs for
      ! control volume field equations such that:
      ! A_m = div(\rho u T) - (1-beta)*T*div(\rho u)

      ! inputs/outputs:
      ! the advection matrix
      type(csr_matrix), intent(inout) :: A_m
      ! the rhs of the control volume field eqn
      type(scalar_field), intent(inout) :: rhs

      ! the field being solved for
      type(scalar_field), intent(inout), target :: tfield
      ! previous time level of the field being solved for
      type(scalar_field), intent(inout) :: oldtfield
      ! a type containing all the tfield options
      type(cv_options_type), intent(in) :: tfield_options
      ! density and previous time level of density associated with the
      ! field (only a real density if solving for
      ! a conservation equation, just constant 1 if AdvectionDiffusion)
      type(scalar_field), intent(inout) :: tdensity, oldtdensity
      ! a type containing all the tdensity options
      type(cv_options_type), intent(in) :: tdensity_options

      ! information about cv faces
      type(cv_faces_type), intent(in) :: cvfaces
      ! shape functions for region and surface
      type(element_type), intent(in) :: x_cvshape, x_cvbdyshape
      type(element_type), intent(in) :: u_cvshape, u_cvbdyshape
      type(element_type), intent(in) :: ug_cvshape, ug_cvbdyshape
      type(element_type), intent(in) :: t_cvshape
      ! bucket full of fields
      type(state_type), dimension(:), intent(inout) :: state
      ! the relative velocity
      type(vector_field), intent(in) :: advu, ug
      ! the coordinates
      type(vector_field), intent(inout) :: x, x_tfield
      ! the cfl number
      type(scalar_field), intent(in) :: cfl_no
      ! logical indicating if the advection matrix should be constructed
      ! or if it exists already from a previous iteration
      logical, intent(in) :: getadvmat
      ! timestep
      real, intent(in) :: dt

      ! mesh sparsity for upwind value matrices
      type(csr_sparsity), intent(in) :: mesh_sparsity

      ! are we using diffusion
      logical, intent(in), optional :: diffusion
      ! logical indicating if the diffusion matrix should be constructed
      ! or if it exists already from a previous iteration
      logical, intent(in), optional :: getdiffmat
      ! the diffusivity tensor
      type(tensor_field), intent(in), optional :: diffusivity
      ! the lumped mass = the mass matrix for the auxilliary diffusion equation
      type(scalar_field), intent(in), optional :: q_lumpedmass
      ! the diffusion matrix
      type(csr_matrix), intent(inout), optional :: D_m
      ! the diffusion rhs
      type(scalar_field), intent(inout), optional :: diff_rhs
      ! sparsity pattern for the gradient transposed operator
      type(csr_sparsity), intent(inout), optional :: grad_m_t_sparsity
      ! shape functions with full body derivatives
      type(element_type), intent(inout), optional :: x_cvshape_full, x_cvbdyshape_full
      type(element_type), intent(inout), optional :: t_cvshape_full, diff_cvshape_full
      type(element_type), intent(inout), optional :: t_cvbdyshape_full, diff_cvbdyshape_full

      ! for back compatibility a reference field with an option path
      ! and appropriate bcs
      type(scalar_field), pointer, optional :: reference_field


      ! local memory:
      ! allocatable memory for coordinates, velocity, normals, determinants, nodes
      ! and the cfl number at the gauss pts and nodes
      real, dimension(:,:), allocatable :: x_ele, x_ele_bdy
      real, dimension(:,:), allocatable :: x_f, u_f, u_bdy_f, ug_f, ug_bdy_f
      real, dimension(:,:), allocatable :: normal, normal_bdy
      real, dimension(:), allocatable :: detwei, detwei_bdy
      real, dimension(:), allocatable :: normgi
      integer, dimension(:), pointer :: nodes, x_nodes, diffusivity_nodes, upwind_nodes
      integer, dimension(:), allocatable :: nodes_bdy, diffusivity_nodes_bdy
      real, dimension(:), allocatable :: cfl_ele

      ! allocatable memory for the values of the field and density at the nodes
      ! and on the boundary and for ghost values outside the boundary
      real, dimension(:), allocatable :: tdensity_ele, oldtdensity_ele, &
                                         tfield_ele, oldtfield_ele
      real, dimension(:), allocatable :: tdensity_ele_bdy, oldtdensity_ele_bdy, &
                                         tfield_ele_bdy, oldtfield_ele_bdy
      real, dimension(:), allocatable :: ghost_tdensity_ele_bdy, ghost_oldtdensity_ele_bdy, &
                                         ghost_tfield_ele_bdy, ghost_gradtfield_ele_bdy, ghost_oldtfield_ele_bdy

      ! some memory used in assembly of the face values
      real :: tfield_theta_val, tdensity_theta_val, tfield_pivot_val
      real :: tfield_face_val, oldtfield_face_val
      real :: tdensity_face_val, oldtdensity_face_val

      ! logical array indicating if a face has already been visited by the opposing node
      logical, dimension(:), allocatable :: notvisited

      ! loop integers
      integer :: ele, sele, iloc, oloc, dloc, face, gi, ggi, dim

      ! upwind value matrices for the fields and densities
      type(csr_matrix)  :: tfield_upwind, &
            oldtfield_upwind, tdensity_upwind, oldtdensity_upwind

      ! incoming or outgoing flow
      real :: udotn, income, udotn_bdy
      logical :: inflow
      ! time and face discretisation
      real :: ptheta, ftheta, beta

      ! the type of the bc if integrating over domain boundaries
      integer, dimension(:), allocatable :: tfield_bc_type, tdensity_bc_type
      ! fields for the bcs over the entire surface mesh
      type(scalar_field) :: tfield_bc, tdensity_bc

      ! back compatible local option path and reference field from the schema
      character(len=OPTION_PATH_LEN) :: l_option_path
      type(scalar_field), pointer :: l_reference_field

      ! local element matrices - allow the assembly of an entire face without multiple calls to csr_pos
      real, dimension(:,:,:), allocatable :: grad_mat_local
      real, dimension(:,:), allocatable :: mat_local, grad_mat_local_bdy, grad_rhs_local_bdy, &
                                           diff_mat_local, diff_mat_local_bdy
      real, dimension(:), allocatable :: mat_local_bdy, rhs_local_bdy, rhs_local, &
                                         div_rhs_local_bdy

      ! a local logical controlling diffusion
      logical :: l_diffusion
      ! the auxilliary gradient matrix (assembled as a divergence confusingly)
      type(block_csr_matrix) :: div_m
      ! the auxilliary gradient equation rhs
      type(vector_field) :: grad_rhs
      ! the diffusivity evaluated at the nodes and the transformed full body gradients
      real, dimension(:,:,:), allocatable :: dt_t, dt_ft, diffusivity_gi, diffusivity_gi_f
      ! a dummy array to potentially store multiple copies of the diffusivity nodes
      integer, dimension(:), allocatable :: diffusivity_lglno, diffusivity_lglno_bdy


      ewrite(2,*) 'assemble_advectiondiffusion_m_cv'

      ! for backward compatibility get the local option path
      ! and field that we need with bcs etc.
      if(present(reference_field)) then
        l_reference_field=>reference_field
        l_option_path = trim(reference_field%option_path)
      else
        l_reference_field=>tfield
        l_option_path = trim(tfield%option_path)
      end if

      if(present(diffusion)) then
        l_diffusion = diffusion
      else
        l_diffusion = .false.
      end if

      if(l_diffusion) then
        select case(tfield_options%diffusionscheme)
        case(CV_DIFFUSION_BASSIREBAY)
          if(.not.present(D_m).or..not.present(diff_rhs).or.&
            .not.present(diffusivity).or..not.present(q_lumpedmass).or.&
            .not.present(grad_m_t_sparsity)) then
            ewrite(-1,*) 'A diffusion matrix and rhs as well as'
            ewrite(-1,*) 'a sparsity, the diffusivity and lumped mass must'
            ewrite(-1,*) 'be supplied to use BassiRebay diffusion!'
            FLAbort("Sorry!")
          end if
        case(CV_DIFFUSION_ELEMENTGRADIENT)
          if(.not.present(D_m).or..not.present(diff_rhs).or.&
            .not.present(diffusivity).or..not.present(x_cvshape_full).or.&
            .not.present(x_cvbdyshape_full).or..not.present(t_cvshape_full).or.&
            .not.present(diff_cvshape_full).or..not.present(t_cvbdyshape_full).or.&
            .not.present(diff_cvbdyshape_full)) then
            ewrite(-1,*) 'A diffusion matrix and rhs as well as'
            ewrite(-1,*) 'the diffusivity and full cv element types must'
            ewrite(-1,*) 'be supplied to use ElementGradient diffusion!'
            FLAbort("Sorry!")
          end if
        end select

        if(present_and_true(getdiffmat)) then
          if(tfield_options%diffusionscheme==CV_DIFFUSION_BASSIREBAY) then
            call allocate(div_m, sparsity=grad_m_t_sparsity, &
                          blocks=(/1, mesh_dim(tfield)/), &
                          name=trim(tfield%name)//"AuxilliaryGradientMatrixTransposed")
            call zero(div_m)
            call allocate(grad_rhs, mesh_dim(tfield), diffusivity%mesh, &
                          name=trim(tfield%name)//"AuxilliaryGradientRHS")
            call zero(grad_rhs)
          end if

          call zero(D_m)
          call zero(diff_rhs)
        end if

      end if

      ! allocate upwind value matrices
      call allocate(tfield_upwind, mesh_sparsity, name="TFieldUpwindValues")
      call allocate(oldtfield_upwind, mesh_sparsity, name="OldTFieldUpwindValues")
      call allocate(tdensity_upwind, mesh_sparsity, name="TDensityUpwindValues")
      call allocate(oldtdensity_upwind, mesh_sparsity, name="OldTDensityUpwindValues")

      ! allocate memory for assembly
      allocate(x_ele(x%dim,ele_loc(x,1)), &
               x_f(x%dim, x_cvshape%ngi), &
               u_f(advu%dim, u_cvshape%ngi), &
               detwei(x_cvshape%ngi), &
               normal(x%dim, x_cvshape%ngi), &
               normgi(x%dim))
      allocate(cfl_ele(ele_loc(cfl_no,1)), &
               tfield_ele(ele_loc(tfield,1)), &
               oldtfield_ele(ele_loc(oldtfield, 1)), &
               tdensity_ele(ele_loc(tdensity, 1)), &
               oldtdensity_ele(ele_loc(oldtdensity,1)))
      allocate(notvisited(x_cvshape%ngi))
      allocate(grad_mat_local(mesh_dim(tfield), tfield%mesh%shape%loc, tfield%mesh%shape%loc), &
               mat_local(tfield%mesh%shape%loc, tfield%mesh%shape%loc), &
               rhs_local(tfield%mesh%shape%loc), &
               diff_mat_local(tfield%mesh%shape%loc, tfield%mesh%shape%loc), &
               dt_t(tfield%mesh%shape%loc, x_cvshape%ngi, mesh_dim(tfield)), &
               diffusivity_gi(mesh_dim(tfield), mesh_dim(tfield), x_cvshape%ngi), &
               diffusivity_lglno(ele_loc(tfield,1)))
      if(move_mesh) then
        allocate(ug_f(ug%dim, ug_cvshape%ngi))
      end if

      ! Clear memory of arrays being designed
      if(getadvmat) call zero(A_m)

      ! does the density field need upwind values?
      if(need_upwind_values(trim(tdensity%option_path))) then

        call find_upwind_values(state, x_tfield, tdensity, tdensity_upwind, &
                                oldtdensity, oldtdensity_upwind &
                                )

      else

        call zero(tdensity_upwind)
        call zero(oldtdensity_upwind)

      end if

      ! does the field need upwind values
      if(need_upwind_values(trim(l_option_path))) then

        call find_upwind_values(state, x_tfield, tfield, tfield_upwind, &
                                oldtfield, oldtfield_upwind, &
                                option_path=trim(l_option_path))

      else

        call zero(tfield_upwind)
        call zero(oldtfield_upwind)

      end if

      ! some temporal discretisation options for clarity
      ptheta = tfield_options%ptheta
      beta = tfield_options%beta

      ! loop over elements
      element_loop: do ele=1, element_count(tfield)
        x_ele=ele_val(x, ele)
        x_f=ele_val_at_quad(x, ele, x_cvshape)
        u_f=ele_val_at_quad(advu, ele, u_cvshape)
        if(move_mesh) ug_f=ele_val_at_quad(ug, ele, ug_cvshape)
        nodes=>ele_nodes(tfield, ele)
        x_nodes=>ele_nodes(x_tfield, ele)
        if((tfield_options%upwind_scheme==CV_UPWINDVALUE_PROJECT_POINT).or.&
           (tfield_options%upwind_scheme==CV_UPWINDVALUE_PROJECT_GRAD)) then
          upwind_nodes=>x_nodes
        else
          upwind_nodes=>nodes
        end if

        ! find determinant and unorientated normal
        call transform_cvsurf_to_physical(x_ele, x_cvshape, &
                                          detwei, normal, cvfaces)

        if(l_diffusion.and.present_and_true(getdiffmat).and.&
            (tfield_options%diffusionscheme==CV_DIFFUSION_BASSIREBAY)) then
          diffusivity_nodes=>ele_nodes(diffusivity, ele)
          ! diffusivity may be on a lower degree mesh than the field... to allow that
          ! without changing the assembly code for each specific case we construct
          ! a mapping to the global nodes that is consistent with the local node
          ! numbering of the parent field.
          ! warning: this is not ideal as it will require more csr_pos's
          ! but its more intended as a proof of concept
          do iloc = 1, size(diffusivity_lglno), size(diffusivity_nodes)
            diffusivity_lglno(iloc:iloc+size(diffusivity_nodes)-1)=diffusivity_nodes
          end do
        else
          diffusivity_nodes=>nodes
          diffusivity_lglno=0
        end if

        if(l_diffusion.and.present_and_true(getdiffmat).and.&
           (tfield_options%diffusionscheme==CV_DIFFUSION_ELEMENTGRADIENT)) then
          call transform_to_physical(X, ele, x_shape=x_cvshape_full, &
                                    shape=t_cvshape_full, dshape=dt_t)
          diffusivity_gi = ele_val_at_quad(diffusivity, ele, diff_cvshape_full)
        else
          dt_t = 0.0
          diffusivity_gi = 0.0
        end if

        cfl_ele = ele_val(cfl_no, ele)

        tfield_ele = ele_val(tfield, ele)
        oldtfield_ele = ele_val(oldtfield, ele)

        tdensity_ele = ele_val(tdensity, ele)
        oldtdensity_ele = ele_val(oldtdensity, ele)

        notvisited=.true.

        grad_mat_local = 0.0
        mat_local = 0.0
        rhs_local = 0.0
        diff_mat_local = 0.0

        ! loop over nodes within this element
        nodal_loop_i: do iloc = 1, tfield%mesh%shape%loc

          ! loop over cv faces internal to this element
          face_loop: do face = 1, cvfaces%faces

            ! is this a face neighbouring iloc?
            if(cvfaces%neiloc(iloc, face) /= 0) then
              oloc = cvfaces%neiloc(iloc, face)

              ! loop over gauss points on face
              quadrature_loop: do gi = 1, cvfaces%shape%ngi

                ! global gauss pt index
                ggi = (face-1)*cvfaces%shape%ngi + gi

                ! have we been here before?
                if(notvisited(ggi)) then
                  notvisited(ggi)=.false.

                  ! correct the orientation of the normal so it points away from iloc
                  normgi=orientate_cvsurf_normgi(node_val(x_tfield, x_nodes(iloc)),x_f(:,ggi),normal(:,ggi))

                  ! calculate u.n
                  if(move_mesh) then
                    udotn=dot_product((u_f(:,ggi)-ug_f(:,ggi)), normgi(:))
                  else
                    udotn=dot_product(u_f(:,ggi), normgi(:))
                  end if
                  inflow = (udotn<=0.0)
                  income = merge(1.0,0.0,inflow)
                  
                  ! calculate the iterated pivot value (so far only does first order upwind)
                  ! which will be subtracted out from the rhs such that with an increasing number
                  ! of iterations the true implicit lhs pivot is cancelled out (if it converges!)
                  tfield_pivot_val = income*tfield_ele(oloc) + (1.-income)*tfield_ele(iloc)

                  ! evaluate the nonlinear face value that will go into the rhs
                  ! this is the value that you choose the discretisation for and
                  ! that will become the dominant term once convergence is achieved
                  call evaluate_face_val(tfield_face_val, oldtfield_face_val, & 
                                         iloc, oloc, ggi, upwind_nodes, &
                                         t_cvshape, &
                                         tfield_ele, oldtfield_ele, &
                                         tfield_upwind, oldtfield_upwind, &
                                         inflow, cfl_ele, &
                                         udotn, &
                                         tfield_options)

                  ! do the same for the density but save some effort if it's just a dummy
                  select case (tdensity%field_type)
                  case(FIELD_TYPE_CONSTANT)

                      tdensity_face_val = tdensity_ele(iloc)
                      oldtdensity_face_val = oldtdensity_ele(iloc)

                  case default

                      call evaluate_face_val(tdensity_face_val, oldtdensity_face_val, &
                                             iloc, oloc, ggi, upwind_nodes, &
                                             t_cvshape,&
                                             tdensity_ele, oldtdensity_ele, &
                                             tdensity_upwind, oldtdensity_upwind, &
                                             inflow, cfl_ele, &
                                             udotn, &
                                             tdensity_options)

                  end select

                  ! perform the time discretisation on the combined tdensity tfield product
                  tfield_theta_val=theta_val(iloc, oloc, &
                                       tfield_face_val, &
                                       oldtfield_face_val, &
                                       tfield_options%theta, dt, udotn, &
                                       x_ele, tfield_options%limit_theta, &
                                       tfield_ele, oldtfield_ele, &
                                       ftheta=ftheta)
                  tdensity_theta_val=theta_val(iloc, oloc, &
                                       tdensity_face_val, &
                                       oldtdensity_face_val, &
                                       tdensity_options%theta, dt, udotn, &
                                       x_ele, tdensity_options%limit_theta, &
                                       tdensity_ele, oldtdensity_ele)

                  if(getadvmat) then
                    mat_local(iloc, oloc) = mat_local(iloc, oloc) &
                                          + ptheta*detwei(ggi)*udotn*income*tdensity_theta_val
                    mat_local(oloc, iloc) = mat_local(oloc, iloc) &
                                          + ptheta*detwei(ggi)*(-udotn)*(1.-income)*tdensity_theta_val
                    mat_local(iloc, iloc) = mat_local(iloc, iloc) &
                                          + ptheta*detwei(ggi)*udotn*(1.0-income)*tdensity_theta_val &
                                          - ftheta*(1.-beta)*detwei(ggi)*udotn*tdensity_theta_val
                    mat_local(oloc, oloc) = mat_local(oloc, oloc) &
                                          + ptheta*detwei(ggi)*(-udotn)*income*tdensity_theta_val &
                                          - ftheta*(1.-beta)*detwei(ggi)*(-udotn)*tdensity_theta_val
                  end if

                  rhs_local(iloc) = rhs_local(iloc) &
                                  + ptheta*udotn*detwei(ggi)*tdensity_theta_val*tfield_pivot_val &
                                  - udotn*detwei(ggi)*tfield_theta_val*tdensity_theta_val &
                                  + (1.-ftheta)*(1.-beta)*detwei(ggi)*udotn*tdensity_theta_val*oldtfield_ele(iloc)
                  rhs_local(oloc) = rhs_local(oloc) &
                                  + ptheta*(-udotn)*detwei(ggi)*tdensity_theta_val*tfield_pivot_val &
                                  - (-udotn)*detwei(ggi)*tfield_theta_val*tdensity_theta_val &
                                  + (1.-ftheta)*(1.-beta)*detwei(ggi)*(-udotn)*tdensity_theta_val*oldtfield_ele(oloc)

                  if(l_diffusion) then

                    if(present_and_true(getdiffmat)) then

                      select case(tfield_options%diffusionscheme)
                      case(CV_DIFFUSION_BASSIREBAY)

                        ! assemble the auxiliary gradient matrix
                        dimension_loop1: do dim = 1, mesh_dim(tfield)

                          grad_mat_local(dim, iloc, iloc) = grad_mat_local(dim, iloc, iloc) &
                                    +0.5*detwei(ggi)*normgi(dim)
!                           grad_mat_local(dim, oloc, iloc) = grad_mat_local(dim, oloc, iloc) &
!                                     +0.5*detwei(ggi)*normgi(dim) ! remember this is a gradient transposed
                          grad_mat_local(dim, iloc, oloc) = grad_mat_local(dim, iloc, oloc) &
                                    +0.5*detwei(ggi)*normgi(dim) ! remember this is a divergence assembly

                          ! notvisited
                          grad_mat_local(dim, oloc, oloc) = grad_mat_local(dim, oloc, oloc) &
                                    -0.5*detwei(ggi)*normgi(dim)
!                           grad_mat_local(dim, iloc, oloc) = grad_mat_local(dim, iloc, oloc) &
!                                     -0.5*detwei(ggi)*normgi(dim) ! remember this is a gradient transposed
                          grad_mat_local(dim, oloc, iloc) = grad_mat_local(dim, oloc, iloc) &
                                    -0.5*detwei(ggi)*normgi(dim) ! remember this is a divergence assembly
!                           ! end of notvisited
                        end do dimension_loop1

                      case(CV_DIFFUSION_ELEMENTGRADIENT)

                        do dloc=1,size(dt_t,1)
                          ! n_i K_{ij} dT/dx_j
                          diff_mat_local(iloc,dloc) = diff_mat_local(iloc,dloc) - &
                            sum(matmul(diffusivity_gi(:,:,ggi), dt_t(dloc, ggi, :))*normgi, 1)*detwei(ggi)

                          ! notvisited
                          diff_mat_local(oloc, dloc) = diff_mat_local(oloc,dloc) - &
                            sum(matmul(diffusivity_gi(:,:,ggi), dt_t(dloc, ggi, :))*(-normgi), 1)*detwei(ggi)
                        end do

                      end select
                    end if

                  end if


                end if ! notvisited
              end do quadrature_loop

            end if ! neiloc
          end do face_loop
        end do nodal_loop_i

        ! if we need the matrix then assemble it now
        if(getadvmat) then
          call addto(A_m, nodes, nodes, mat_local)
        end if

        if(present_and_true(getdiffmat)) then
          if(l_diffusion) then
            select case(tfield_options%diffusionscheme)
            case(CV_DIFFUSION_BASSIREBAY)

              call addto(div_m, nodes, diffusivity_lglno, spread(grad_mat_local, 1, 1))

            case(CV_DIFFUSION_ELEMENTGRADIENT)

              call addto(D_m, nodes, nodes, diff_mat_local)

            end select
          end if
        end if

        ! assemble the rhs
        call addto(rhs, nodes, rhs_local)

      end do element_loop

      ! allocate memory for assembly
      allocate(x_ele_bdy(x%dim,face_loc(x,1)), &
              detwei_bdy(x_cvbdyshape%ngi), &
              normal_bdy(x%dim, x_cvbdyshape%ngi), &
              u_bdy_f(advu%dim, u_cvbdyshape%ngi), &
              tdensity_ele_bdy(face_loc(tdensity,1)), &
              oldtdensity_ele_bdy(face_loc(oldtdensity,1)), &
              tfield_ele_bdy(face_loc(tfield,1)), &
              oldtfield_ele_bdy(face_loc(oldtfield,1)), &
              ghost_tdensity_ele_bdy(face_loc(tdensity,1)), &
              ghost_oldtdensity_ele_bdy(face_loc(oldtdensity,1)), &
              ghost_tfield_ele_bdy(face_loc(tfield,1)), &
              ghost_gradtfield_ele_bdy(face_loc(tfield,1)), &
              ghost_oldtfield_ele_bdy(face_loc(oldtfield,1)))
      allocate(tfield_bc_type(surface_element_count(tfield)), &
              tdensity_bc_type(surface_element_count(tdensity)), &
              nodes_bdy(face_loc(tfield,1)))
      allocate(grad_mat_local_bdy(mesh_dim(tfield), tfield%mesh%faces%shape%loc), &
               grad_rhs_local_bdy(mesh_dim(tfield), tfield%mesh%faces%shape%loc), &
               div_rhs_local_bdy(tfield%mesh%faces%shape%loc), &
               mat_local_bdy(tfield%mesh%faces%shape%loc), &
               rhs_local_bdy(tfield%mesh%faces%shape%loc), &
               dt_ft(tfield%mesh%faces%shape%loc, x_cvbdyshape%ngi, mesh_dim(tfield)), &
               diffusivity_gi_f(mesh_dim(tfield), mesh_dim(tfield), x_cvbdyshape%ngi), &
               diff_mat_local_bdy(tfield%mesh%faces%shape%loc, tfield%mesh%faces%shape%loc), &
               diffusivity_lglno_bdy(face_loc(tfield,1)))
      if(l_diffusion) then
        allocate(diffusivity_nodes_bdy(face_loc(diffusivity,1)))
      else
        allocate(diffusivity_nodes_bdy(face_loc(tfield,1)))
      end if
      if(move_mesh) then
        allocate(ug_bdy_f(ug%dim, ug_cvbdyshape%ngi))
      end if

      ! get the fields over the surface containing the bcs
      call get_entire_boundary_condition(l_reference_field, (/ &
        "weakdirichlet", &
        "neumann      ", &
        "periodic     ", &
        "zero_flux    "/), tfield_bc, tfield_bc_type)
      call get_entire_boundary_condition(tdensity, (/"weakdirichlet"/), tdensity_bc, tdensity_bc_type)

      ! loop over the surface elements
      surface_element_loop: do sele = 1, surface_element_count(tfield)
        
        if((tfield_bc_type(sele)==3)) cycle

        ele = face_ele(x, sele)
        x_ele = ele_val(x, ele)
        x_ele_bdy = face_val(x, sele)
        nodes_bdy=face_global_nodes(tfield, sele)

        ! calculate the determinant and orientated normal
        call transform_cvsurf_facet_to_physical(x_ele, x_ele_bdy, &
                              x_cvbdyshape, normal_bdy, detwei_bdy)

        if(l_diffusion.and.present_and_true(getdiffmat).and.&
            (tfield_options%diffusionscheme==CV_DIFFUSION_BASSIREBAY)) then
          diffusivity_nodes_bdy=face_global_nodes(diffusivity,sele)
          ! diffusivity may be on a lower degree mesh than the field... to allow that
          ! without changing the assembly code for each specific case we construct
          ! a mapping to the global nodes that is consistent with the local node
          ! numbering of the parent field.
          ! warning: this is not ideal as it will require more csr_pos's
          ! but its more intended as a proof of concept
          do iloc = 1, size(diffusivity_lglno_bdy), size(diffusivity_nodes_bdy)
            diffusivity_lglno_bdy(iloc:iloc+size(diffusivity_nodes_bdy)-1)=diffusivity_nodes_bdy
          end do
        else
          diffusivity_nodes_bdy=0
          diffusivity_lglno_bdy=0
        end if

        if(l_diffusion.and.present_and_true(getdiffmat).and.&
           (tfield_options%diffusionscheme==CV_DIFFUSION_ELEMENTGRADIENT)) then
!           call transform_to_physical(x_ele_bdy, x_cvbdyshape_full, &
!                                     m=t_cvbdyshape_full, dm_t=dt_ft)
          dt_ft = 0.0 ! at the moment its not possible to get the full gradient
                      ! so until this is fixed we're just going to have to assume 
                      ! zero neumann on outflow boundaries
          diffusivity_gi_f = face_val_at_quad(diffusivity, sele, diff_cvbdyshape_full)
        else
          dt_ft = 0.0
          diffusivity_gi_f = 0.0
        end if

        u_bdy_f=face_val_at_quad(advu, sele, u_cvbdyshape)
        if(move_mesh) ug_bdy_f=face_val_at_quad(ug, sele, ug_cvbdyshape)

        ! deal with bcs for tfield
        if(tfield_bc_type(sele)==1) then
          ghost_tfield_ele_bdy=ele_val(tfield_bc, sele)
        else
          ghost_tfield_ele_bdy=face_val(tfield, sele)
        end if

        ghost_gradtfield_ele_bdy = ele_val(tfield_bc, sele)

        if(tfield_bc_type(sele)==1) then
          ghost_oldtfield_ele_bdy=ele_val(tfield_bc, sele) ! not considering time varying bcs yet
        else
          ghost_oldtfield_ele_bdy=face_val(oldtfield, sele)
        end if

        tfield_ele_bdy=face_val(tfield, sele)
        oldtfield_ele_bdy=face_val(oldtfield, sele)

        ! deal with bcs for tdensity
        if(tdensity_bc_type(sele)==1) then
          ghost_tdensity_ele_bdy=ele_val(tdensity_bc, sele)
        else
          ghost_tdensity_ele_bdy=face_val(tdensity, sele)
        end if

        if(tdensity_bc_type(sele)==1) then
          ghost_oldtdensity_ele_bdy=ele_val(tdensity_bc, sele) ! not considering time varying bcs yet
        else
          ghost_oldtdensity_ele_bdy=face_val(oldtdensity, sele)
        end if

        tdensity_ele_bdy=face_val(tdensity, sele)
        oldtdensity_ele_bdy=face_val(oldtdensity, sele)

        grad_mat_local_bdy = 0.0
        grad_rhs_local_bdy = 0.0
        div_rhs_local_bdy = 0.0
        mat_local_bdy = 0.0
        rhs_local_bdy = 0.0
        diff_mat_local_bdy = 0.0

        ! loop over the nodes on this surface element
        surface_nodal_loop_i: do iloc = 1, tfield%mesh%faces%shape%loc

          ! loop over the faces in this surface element
          surface_face_loop: do face = 1, cvfaces%sfaces

            ! is this face a neighbour of iloc?
            if(cvfaces%sneiloc(iloc,face)/=0) then

              ! loop over the gauss pts on this face
              surface_quadrature_loop: do gi = 1, cvfaces%shape%ngi

                ! global gauss point index
                ggi = (face-1)*cvfaces%shape%ngi + gi

                ! u.n
                if(move_mesh) then
                  udotn_bdy = dot_product((u_bdy_f(:,ggi)-ug_bdy_f(:,ggi)), normal_bdy(:,ggi))
                else
                  udotn_bdy = dot_product(u_bdy_f(:,ggi), normal_bdy(:,ggi))
                end if
                
                if((tfield_bc_type(sele)==4)) then
                  ! zero_flux
                  udotn=0.0
                else
                  udotn=udotn_bdy
                end if

                if(udotn>0) then
                  income=0.0 ! flow leaving the domain
                else
                  income=1.0 ! flow entering the domain
                end if

                ! as we're on the boundary it's not possible to use "high order" methods so just
                ! default to the pivotted solution method (first order upwinding)
                ! if the flow is incoming then use the bc ghost values
                ! if the flow is outgoing then use the surface nodes value

                ! for tfield
                tfield_face_val = income*ghost_tfield_ele_bdy(iloc) + (1.-income)*tfield_ele_bdy(iloc)
                oldtfield_face_val = income*ghost_oldtfield_ele_bdy(iloc) + (1.-income)*oldtfield_ele_bdy(iloc)

                ! for tdensity
                tdensity_face_val = income*ghost_tdensity_ele_bdy(iloc) + (1.-income)*tdensity_ele_bdy(iloc)
                oldtdensity_face_val = income*ghost_oldtdensity_ele_bdy(iloc) + (1.-income)*oldtdensity_ele_bdy(iloc)

                tdensity_theta_val = tdensity_options%theta*tdensity_face_val + (1.-tdensity_options%theta)*oldtdensity_face_val

                if(getadvmat) then
                  ! if iloc is the donor we can do this implicitly
                  mat_local_bdy(iloc) = mat_local_bdy(iloc) &
                                  + ptheta*detwei_bdy(ggi)*udotn*(1.-income)*tdensity_theta_val &  
                                  - ptheta*(1.-beta)*detwei_bdy(ggi)*udotn_bdy*tdensity_theta_val
                end if

                ! but we can't if it's the downwind
                rhs_local_bdy(iloc) = rhs_local_bdy(iloc) &
                              - ptheta*udotn*detwei_bdy(ggi)*income*tdensity_theta_val*ghost_tfield_ele_bdy(iloc) & 
                              - (1.-ptheta)*udotn*detwei_bdy(ggi)*tdensity_theta_val*oldtfield_face_val &
                              + (1.-ptheta)*(1.-beta)*udotn_bdy*detwei_bdy(ggi)*tdensity_theta_val*oldtfield_ele_bdy(iloc)

                if(l_diffusion.and.present_and_true(getdiffmat)) then

                  select case(tfield_options%diffusionscheme)
                  case(CV_DIFFUSION_BASSIREBAY)

                    if(tfield_bc_type(sele)==1) then
                      ! assemble grad_rhs

                      grad_rhs_local_bdy(:, iloc) = grad_rhs_local_bdy(:,iloc) &
                                  -detwei_bdy(ggi)*normal_bdy(:,ggi)*ghost_tfield_ele_bdy(iloc)

                      ! when assembling a divergence operator you need this:
                      ! (but not when its a gradient transposed operator)
                      dimension_loop2: do dim = 1, mesh_dim(tfield)
                        ! assemble matrix
                        grad_mat_local_bdy(dim, iloc) = grad_mat_local_bdy(dim, iloc) &
                                          +detwei_bdy(ggi)*normal_bdy(dim,ggi)
                      end do dimension_loop2

                    else

                      if(tfield_bc_type(sele)==2) then

                        ! assemble div_rhs
                        div_rhs_local_bdy(iloc) = div_rhs_local_bdy(iloc) &
                                    -detwei_bdy(ggi)*ghost_gradtfield_ele_bdy(iloc)

                      end if


!                       ! when assembling a gradient transposed operator you need this:
!                       ! (but not when its a divergence operator)
!                       dimension_loop2: do dim = 1, mesh_dim(tfield)
!                         ! assemble matrix
!                         grad_mat_local_bdy(dim, iloc) = grad_mat_local_bdy(dim, iloc) &
!                                           -detwei_bdy(ggi)*normal_bdy(dim,ggi)
!                       end do dimension_loop2

                    end if

                  case(CV_DIFFUSION_ELEMENTGRADIENT)

                    if(tfield_bc_type(sele)==2) then

                      div_rhs_local_bdy(iloc) = div_rhs_local_bdy(iloc) &
                                  -detwei_bdy(ggi)*ghost_gradtfield_ele_bdy(iloc)

                    else

!                 because transform to physical doesn't give the full gradient at a face
!                 yet this can't be done so we're going to have to assume zero neumann
!                 at outflow faces
!                       do dloc= 1,tfield%mesh%faces%shape%loc
! 
!                         ! n_i K_{ij} dT/dx_j
!                         diff_mat_local_bdy(iloc, dloc) = diff_mat_local_bdy(iloc,dloc) + &
!                           sum(matmul(diffusivity_gi_f(:,:,ggi), dt_ft(dloc, ggi, :))*normal_bdy(:,ggi), 1)&
!                           *detwei_bdy(ggi)
! 
!                       end do

                    end if

                  end select

                end if

              end do surface_quadrature_loop

            end if ! sneiloc

          end do surface_face_loop

        end do surface_nodal_loop_i

        ! assemble matrix
        if(getadvmat) then
          call addto_diag(A_m, nodes_bdy, mat_local_bdy)
        end if

        if(present_and_true(getdiffmat)) then
          if(l_diffusion) then
            select case(tfield_options%diffusionscheme)
            case(CV_DIFFUSION_BASSIREBAY)

              do dim = 1, mesh_dim(tfield)
                do iloc = 1, size(grad_mat_local_bdy,2)
                  call addto(div_m, 1, dim, nodes_bdy(iloc), diffusivity_lglno_bdy(iloc), &
                             grad_mat_local_bdy(dim,iloc))
                end do
              end do
              call addto(grad_rhs, diffusivity_lglno_bdy, grad_rhs_local_bdy)

              call addto(diff_rhs, nodes_bdy, div_rhs_local_bdy)

            case(CV_DIFFUSION_ELEMENTGRADIENT)

              if(tfield_bc_type(sele)==1) then

!               assume zero neumann for the moment
!                 call addto(diff_rhs, nodes_bdy, -matmul(diff_mat_local_bdy, ghost_gradtfield_ele_bdy))

              elseif(tfield_bc_type(sele)==2) then

                call addto(diff_rhs, nodes_bdy, div_rhs_local_bdy)

              else

!               assume zero neumann for the moment
!                 call addto(D_m, nodes_bdy, nodes_bdy, diff_mat_local_bdy)

              end if

            end select
          end if
        end if

        ! assemble rhs
        call addto(rhs, nodes_bdy, rhs_local_bdy)

      end do surface_element_loop

      if(present_and_true(getdiffmat).and.l_diffusion.and.&
        (tfield_options%diffusionscheme==CV_DIFFUSION_BASSIREBAY)) then

        call assemble_bassirebay_diffusion_m_cv(D_m, diff_rhs, &
                                     div_m, grad_rhs, &
                                     diffusivity, q_lumpedmass)

      end if

      deallocate(x_ele_bdy, detwei_bdy, normal_bdy, u_bdy_f)
      if(move_mesh) deallocate(ug_bdy_f)
      deallocate(nodes_bdy)
      deallocate(tdensity_ele_bdy, oldtdensity_ele_bdy, tfield_ele_bdy, oldtfield_ele_bdy)
      deallocate(ghost_tdensity_ele_bdy, ghost_oldtdensity_ele_bdy, &
                  ghost_tfield_ele_bdy, ghost_gradtfield_ele_bdy, ghost_oldtfield_ele_bdy)

      deallocate(tfield_bc_type, tdensity_bc_type)
      call deallocate(tfield_bc)
      call deallocate(tdensity_bc)

      call deallocate(tdensity_upwind)
      call deallocate(oldtdensity_upwind)
      call deallocate(tfield_upwind)
      call deallocate(oldtfield_upwind)

      deallocate(x_ele, x_f, detwei, normal, normgi, u_f)
      deallocate(cfl_ele, tfield_ele, oldtfield_ele, tdensity_ele, oldtdensity_ele)
      deallocate(notvisited)

      if(present_and_true(getdiffmat).and.l_diffusion.and.&
        (tfield_options%diffusionscheme==CV_DIFFUSION_BASSIREBAY)) then
        call deallocate(div_m)
        call deallocate(grad_rhs)
      end if

    end subroutine assemble_advectiondiffusion_m_cv

    subroutine assemble_coupled_advection_m_cv(A_m, rhs, &
                                       tfield, oldtfield, tfield_options, &
                                       tdensity, oldtdensity, tdensity_options, &
                                       cvfaces, x_cvshape, x_cvbdyshape, &
                                       u_cvshape, u_cvbdyshape, t_cvshape, &
                                       state, advu, x, x_tfield, cfl_no, getmat, dt, &
                                       mesh_sparsity)

      !!< This subroutine assembles the advection matrix and rhs for
      !!< control volume field equations such that:
      !!< A_m = div(\rho u T) - (1-beta)*T*div(\rho u)
      !!< with the added restrictions between different material's face values

      ! inputs/outputs:
      ! the advection matrix
      type(csr_matrix), dimension(:), intent(inout) :: A_m
      ! the rhs of the control volume field eqn
      type(scalar_field), dimension(:), intent(inout) :: rhs

      ! the fields being solved for
      type(scalar_field_pointer), dimension(:), intent(inout) :: tfield
      ! previous time level of the fields being solved for
      type(scalar_field_pointer), dimension(:), intent(inout) :: oldtfield
      ! a type containing all the tfield options
      type(cv_options_type), dimension(:), intent(in) :: tfield_options
      ! density and previous time level of density associated with the
      ! field (only a real density if solving for
      ! a conservation equation, just constant 1 if AdvectionDiffusion)
      type(scalar_field_pointer), dimension(:), intent(inout) :: tdensity, oldtdensity
      ! a type containing all the tdensity options
      type(cv_options_type), dimension(:), intent(in) :: tdensity_options

      ! information about cv faces
      type(cv_faces_type), intent(in) :: cvfaces
      ! shape functions for region and surface
      type(element_type), intent(in) :: x_cvshape, x_cvbdyshape
      type(element_type), intent(in) :: u_cvshape, u_cvbdyshape
      type(element_type), intent(in) :: t_cvshape
      ! bucket full of fields
      type(state_type), dimension(:), intent(inout) :: state
      ! the advection velocity
      type(vector_field), intent(in) :: advu
      ! the coordinates
      type(vector_field), intent(inout) :: x, x_tfield
      ! the cfl number
      type(scalar_field), intent(in) :: cfl_no
      ! logical indicating if the matrix should be constructed
      ! or if it exists already from a previous iteration
      logical, dimension(:), intent(in) :: getmat
      ! timestep
      real, intent(in) :: dt

      ! mesh sparsity for upwind value matrices
      type(csr_sparsity), intent(in) :: mesh_sparsity

      ! local memory:
      ! allocatable memory for coordinates, velocity, normals, determinants, nodes
      ! and the cfl number at the gauss pts and nodes
      real, dimension(:,:), allocatable :: x_ele, x_ele_bdy
      real, dimension(:,:), allocatable :: x_f, u_f, u_bdy_f
      real, dimension(:,:), allocatable :: normal, normal_bdy
      real, dimension(:), allocatable :: detwei, detwei_bdy
      real, dimension(:), allocatable :: normgi
      integer, dimension(:), pointer :: nodes, x_nodes, upwind_nodes
      integer, dimension(:), allocatable :: nodes_bdy
      real, dimension(:), allocatable :: cfl_ele

      ! allocatable memory for the values of the field and density at the nodes
      ! and on the boundary and for ghost values outside the boundary
      real, dimension(:,:), allocatable :: tdensity_ele, oldtdensity_ele, &
                                           tfield_ele, oldtfield_ele, &
                                           sum_tfield_ele, sum_oldtfield_ele
      real, dimension(:,:), allocatable :: tdensity_ele_bdy, oldtdensity_ele_bdy, &
                                           tfield_ele_bdy, oldtfield_ele_bdy
      real, dimension(:,:), allocatable :: ghost_tdensity_ele_bdy, ghost_oldtdensity_ele_bdy, &
                                           ghost_tfield_ele_bdy, ghost_oldtfield_ele_bdy

      ! some memory used in assembly of the face values
      real :: tfield_theta_val, tfield_pivot_val, tdensity_theta_val
      real, dimension(size(tfield)) :: tfield_face_val, oldtfield_face_val
      real :: tdensity_face_val, oldtdensity_face_val

      ! logical array indicating if a face has already been visited by the opposing node
      logical, dimension(:), allocatable :: notvisited

      ! loop integers
      integer :: ele, sele, iloc, oloc, face, gi, ggi

      ! upwind value matrices for the fields and densities
      type(csr_matrix), dimension(size(tfield)) :: tfield_upwind, &
            oldtfield_upwind, tdensity_upwind, oldtdensity_upwind

      ! incoming or outgoing flow
      real :: udotn, income, udotn_bdy
      logical :: inflow
      ! time and face discretisation
      real, dimension(size(tfield)) :: ptheta, beta
      real :: ftheta

      ! the type of the bc if integrating over domain boundaries
      integer, dimension(:,:), allocatable :: tfield_bc_type, tdensity_bc_type
      ! fields for the bcs over the entire surface mesh
      type(scalar_field), dimension(size(tfield)) :: tfield_bc, tdensity_bc

      integer :: f, f2, nfields, upwind_pos

      ewrite(2,*) 'in assemble_coupled_advection_m_cv'

      nfields = size(tfield)
      upwind_pos = 0

      ! allocate upwind value matrices
      do f = 1, nfields
        call allocate(tfield_upwind(f), mesh_sparsity, name=int2str(f)//"TFieldUpwindValues")
        call allocate(oldtfield_upwind(f), mesh_sparsity, name=int2str(f)//"OldTFieldUpwindValues")
        call allocate(tdensity_upwind(f), mesh_sparsity, name=int2str(f)//"TDensityUpwindValues")
        call allocate(oldtdensity_upwind(f), mesh_sparsity, name=int2str(f)//"OldTDensityUpwindValues")
      end do

      ! allocate memory for assembly
      allocate(x_ele(x%dim,ele_loc(x,1)), &
               x_f(x%dim, x_cvshape%ngi), &
               u_f(advu%dim, u_cvshape%ngi), &
               detwei(x_cvshape%ngi), &
               normal(x%dim, x_cvshape%ngi), &
               normgi(x%dim))
      allocate(cfl_ele(ele_loc(cfl_no,1)), &
               tfield_ele(nfields, ele_loc(tfield(1)%ptr,1)), &
               oldtfield_ele(nfields, ele_loc(oldtfield(1)%ptr, 1)), &
               sum_tfield_ele(nfields, ele_loc(tfield(1)%ptr,1)), &
               sum_oldtfield_ele(nfields, ele_loc(oldtfield(1)%ptr, 1)), &
               tdensity_ele(nfields, ele_loc(tdensity(1)%ptr, 1)), &
               oldtdensity_ele(nfields, ele_loc(oldtdensity(1)%ptr,1)))
      allocate(notvisited(x_cvshape%ngi))

      ! Clear memory of arrays being designed
      do f = 1, nfields
        if(getmat(f)) then
          call zero(A_m(f))
        end if
      end do

      ! does the density field need upwind values?
      do f = 1, nfields
        if(need_upwind_values(trim(tdensity(f)%ptr%option_path))) then

          call find_upwind_values(state, x_tfield, tdensity(f)%ptr, tdensity_upwind(f), &
                                  oldtdensity(f)%ptr, oldtdensity_upwind(f) &
                                  )

        else

          call zero(tdensity_upwind(f))
          call zero(oldtdensity_upwind(f))

        end if

        ! does the field need upwind values
        if(need_upwind_values(trim(tfield(f)%ptr%option_path))) then

          call find_upwind_values(state, x_tfield, tfield(f)%ptr, tfield_upwind(f), &
                                  oldtfield(f)%ptr, oldtfield_upwind(f))

        else

          call zero(tfield_upwind(f))
          call zero(oldtfield_upwind(f))

        end if

        ! some temporal discretisation options for clarity
        ptheta(f) = tfield_options(f)%ptheta
        beta(f) = tfield_options(f)%beta
        
      end do
      
      call couple_upwind_values(tfield_upwind, oldtfield_upwind, tfield_options)
      
      ! loop over elements
      element_loop: do ele=1, element_count(tfield(1)%ptr)
        x_ele=ele_val(x, ele)
        x_f=ele_val_at_quad(x, ele, x_cvshape)
        u_f=ele_val_at_quad(advu, ele, u_cvshape)
        nodes=>ele_nodes(tfield(1)%ptr, ele)
        x_nodes=>ele_nodes(x_tfield, ele)
        if((tfield_options(1)%upwind_scheme==CV_UPWINDVALUE_PROJECT_POINT).or.&
           (tfield_options(1)%upwind_scheme==CV_UPWINDVALUE_PROJECT_GRAD)) then
          upwind_nodes=>x_nodes
        else
          upwind_nodes=>nodes
        end if

        ! find determinant and unorientated normal
        call transform_cvsurf_to_physical(x_ele, x_cvshape, &
                                          detwei, normal, cvfaces)

        cfl_ele = ele_val(cfl_no, ele)

        sum_tfield_ele = 0.0
        sum_oldtfield_ele = 0.0
        do f = 1, nfields
          tfield_ele(f,:) = ele_val(tfield(f)%ptr, ele)
          oldtfield_ele(f,:) = ele_val(oldtfield(f)%ptr, ele)
          
          do f2 = f, nfields
            sum_tfield_ele(f2,:) = sum_tfield_ele(f2,:) + tfield_ele(f,:)
            sum_oldtfield_ele(f2,:) = sum_oldtfield_ele(f2,:) + oldtfield_ele(f,:)
          end do

          tdensity_ele(f,:) = ele_val(tdensity(f)%ptr, ele)
          oldtdensity_ele(f,:) = ele_val(oldtdensity(f)%ptr, ele)
        end do

        notvisited=.true.

        ! loop over nodes within this element
        nodal_loop_i: do iloc = 1, tfield(1)%ptr%mesh%shape%loc

          ! loop over cv faces internal to this element
          face_loop: do face = 1, cvfaces%faces
          
            ! is this a face neighbouring iloc?
            if(cvfaces%neiloc(iloc, face) /= 0) then
              oloc = cvfaces%neiloc(iloc, face)

              ! loop over gauss points on face
              quadrature_loop: do gi = 1, cvfaces%shape%ngi

                ! global gauss pt index
                ggi = (face-1)*cvfaces%shape%ngi + gi

                ! have we been here before?
                if(notvisited(ggi)) then
                  notvisited(ggi)=.false.

                  ! correct the orientation of the normal so it points away from iloc
                  normgi=orientate_cvsurf_normgi(node_val(x_tfield, x_nodes(iloc)),x_f(:,ggi),normal(:,ggi))

                  ! calculate u.n
                  udotn=dot_product(u_f(:,ggi), normgi(:))

                  inflow = (udotn<=0.0)

                  income = merge(1.0,0.0,inflow)

                  field_loop: do f = 1, nfields
                    ! calculate the iterated pivot value (so far only does first order upwind)
                    ! which will be subtracted out from the rhs such that with an increasing number
                    ! of iterations the true implicit lhs pivot is cancelled out (if it converges!)
                    tfield_pivot_val = income*tfield_ele(f, oloc) + (1.-income)*tfield_ele(f, iloc)

                    ! evaluate the nonlinear face value that will go into the rhs
                    ! this is the value that you choose the discretisation for and
                    ! that will become the dominant term once convergence is achieved

                    ! do the same for the density but save some effort if it's just a dummy
                    select case (tdensity(f)%ptr%field_type)
                    case(FIELD_TYPE_CONSTANT)

                        tdensity_face_val = tdensity_ele(f,iloc)
                        oldtdensity_face_val = oldtdensity_ele(f,iloc)

                    case default

                        call evaluate_face_val(tdensity_face_val, oldtdensity_face_val, &
                                              iloc, oloc, ggi, upwind_nodes, &
                                              t_cvshape,&
                                              tdensity_ele(f,:), oldtdensity_ele(f,:), &
                                              tdensity_upwind(f), oldtdensity_upwind(f), &
                                              inflow, cfl_ele, &
                                              udotn, &
                                              tdensity_options(f), save_pos = upwind_pos)

                    end select

                    call evaluate_face_val(tfield_face_val(f), oldtfield_face_val(f), & 
                                          iloc, oloc, ggi, upwind_nodes, &
                                          t_cvshape, &
                                          tfield_ele(f,:), oldtfield_ele(f,:), &
                                          tfield_upwind(f), oldtfield_upwind(f), &
                                          inflow, cfl_ele, &
                                          udotn, &
                                          tfield_options(f), save_pos=upwind_pos)
                    
                    if(f>1) then
                    
                      call couple_face_value(tfield_face_val(f), oldtfield_face_val(f), &
                                             sum(tfield_face_val(1:f-1)), sum(oldtfield_face_val(1:f-1)), &
                                             tfield_ele(f,:), oldtfield_ele(f,:), &
                                             sum_tfield_ele(f-1,:), sum_oldtfield_ele(f-1,:), &
                                             tfield_upwind(f), oldtfield_upwind(f), &
                                             inflow, iloc, oloc, upwind_nodes, cfl_ele, &
                                             tfield_options(f), save_pos=upwind_pos)

                    end if

                    ! perform the time discretisation on the combined tdensity tfield product
                    tfield_theta_val=theta_val(iloc, oloc, &
                                        tfield_face_val(f), &
                                        oldtfield_face_val(f), &
                                        tfield_options(f)%theta, dt, udotn, &
                                        x_ele, tfield_options(f)%limit_theta, &
                                        tfield_ele(f,:), oldtfield_ele(f,:), &
                                        ftheta=ftheta)

                    tdensity_theta_val=theta_val(iloc, oloc, &
                                        tdensity_face_val, &
                                        oldtdensity_face_val, &
                                        tdensity_options(f)%theta, dt, udotn, &
                                        x_ele, tdensity_options(f)%limit_theta, &
                                        tdensity_ele(f,:), oldtdensity_ele(f,:))

                    ! if we need the matrix then assemble it now
                    if(getmat(f)) then
                      call addto(A_m(f), nodes(iloc), nodes(oloc), &
                                ptheta(f)*detwei(ggi)*udotn*income*tdensity_theta_val)
                      call addto(A_m(f), nodes(oloc), nodes(iloc), &
                                ptheta(f)*detwei(ggi)*(-udotn)*(1.-income)*tdensity_theta_val) ! notvisited

                      call addto_diag(A_m(f), nodes(iloc), &
                                ptheta(f)*detwei(ggi)*udotn*(1.0-income)*tdensity_theta_val &
                                -ftheta*(1.-beta(f))*detwei(ggi)*udotn*tdensity_theta_val)
                      call addto_diag(A_m(f), nodes(oloc), &
                                ptheta(f)*detwei(ggi)*(-udotn)*income*tdensity_theta_val &
                                -ftheta*(1.-beta(f))*detwei(ggi)*(-udotn)*tdensity_theta_val) ! notvisited

                    end if

                    ! assemble the rhs
                    call addto(rhs(f), nodes(iloc), &
                                  ptheta(f)*udotn*detwei(ggi)*tdensity_theta_val*tfield_pivot_val &
                                - udotn*detwei(ggi)*tfield_theta_val*tdensity_theta_val &
                                + (1.-ftheta)*(1.-beta(f))*detwei(ggi)*udotn*tdensity_theta_val*oldtfield_ele(f,iloc))
                    call addto(rhs(f), nodes(oloc), &
                                  ptheta(f)*(-udotn)*detwei(ggi)*tdensity_theta_val*tfield_pivot_val &
                                - (-udotn)*detwei(ggi)*tfield_theta_val*tdensity_theta_val &
                                + (1.-ftheta)*(1.-beta(f))*detwei(ggi)*(-udotn)*tdensity_theta_val*oldtfield_ele(f,oloc)) ! notvisited

                  end do field_loop
                  
                end if ! notvisited
              end do quadrature_loop
            end if ! neiloc
          end do face_loop
        end do nodal_loop_i
      end do element_loop
      
      ! allocate memory for assembly
      allocate(x_ele_bdy(x%dim,face_loc(x,1)), &
              detwei_bdy(x_cvbdyshape%ngi), &
              normal_bdy(x%dim, x_cvbdyshape%ngi), &
              u_bdy_f(advu%dim, u_cvbdyshape%ngi), &
              tdensity_ele_bdy(nfields,face_loc(tdensity(1)%ptr,1)), &
              oldtdensity_ele_bdy(nfields,face_loc(oldtdensity(1)%ptr,1)), &
              tfield_ele_bdy(nfields,face_loc(tfield(1)%ptr,1)), &
              oldtfield_ele_bdy(nfields,face_loc(oldtfield(1)%ptr,1)), &
              ghost_tdensity_ele_bdy(nfields,face_loc(tdensity(1)%ptr,1)), &
              ghost_oldtdensity_ele_bdy(nfields,face_loc(oldtdensity(1)%ptr,1)), &
              ghost_tfield_ele_bdy(nfields,face_loc(tfield(1)%ptr,1)), &
              ghost_oldtfield_ele_bdy(nfields,face_loc(oldtfield(1)%ptr,1)))
      allocate(tfield_bc_type(nfields, surface_element_count(tfield(1)%ptr)), &
              tdensity_bc_type(nfields, surface_element_count(tdensity(1)%ptr)), &
              nodes_bdy(face_loc(tfield(1)%ptr,1)))

      do f = 1, nfields
        ! get the fields over the surface containing the bcs
        call get_entire_boundary_condition(tfield(f)%ptr, (/ &
          "weakdirichlet", &
          "periodic     ", &
          "zero_flux    "/), tfield_bc(f), tfield_bc_type(f,:))
        call get_entire_boundary_condition(tdensity(f)%ptr, (/"weakdirichlet"/), tdensity_bc(f), tdensity_bc_type(f,:))
      end do
      
      ! loop over the surface elements
      surface_element_loop: do sele = 1, surface_element_count(tfield(1)%ptr)

        if(any(tfield_bc_type(:,sele)==2)) cycle
        
        ele = face_ele(x, sele)
        x_ele = ele_val(x, ele)
        x_ele_bdy = face_val(x, sele)
        nodes_bdy=face_global_nodes(tfield(1)%ptr, sele)

        ! calculate the determinant and orientated normal
        call transform_cvsurf_facet_to_physical(x_ele, x_ele_bdy, &
                              x_cvbdyshape, normal_bdy, detwei_bdy)

        u_bdy_f=face_val_at_quad(advu, sele, u_cvbdyshape)

        do f = 1, nfields
          ! deal with bcs for tfield
          if(tfield_bc_type(f,sele)==1) then
            ghost_tfield_ele_bdy(f,:)=ele_val(tfield_bc(f), sele)
          else
            ghost_tfield_ele_bdy(f,:)=face_val(tfield(f)%ptr, sele)
          end if

          if(tfield_bc_type(f,sele)==1) then
            ghost_oldtfield_ele_bdy(f,:)=ele_val(tfield_bc(f), sele) ! not considering time varying bcs yet
          else
            ghost_oldtfield_ele_bdy(f,:)=face_val(oldtfield(f)%ptr, sele)
          end if

          tfield_ele_bdy(f,:)=face_val(tfield(f)%ptr, sele)
          oldtfield_ele_bdy(f,:)=face_val(oldtfield(f)%ptr, sele)

          ! deal with bcs for tdensity
          if(tdensity_bc_type(f,sele)==1) then
            ghost_tdensity_ele_bdy(f,:)=ele_val(tdensity_bc(f), sele)
          else
            ghost_tdensity_ele_bdy(f,:)=face_val(tdensity(f)%ptr, sele)
          end if

          if(tdensity_bc_type(f,sele)==1) then
            ghost_oldtdensity_ele_bdy(f,:)=ele_val(tdensity_bc(f), sele) ! not considering time varying bcs yet
          else
            ghost_oldtdensity_ele_bdy(f,:)=face_val(oldtdensity(f)%ptr, sele)
          end if

          tdensity_ele_bdy(f,:)=face_val(tdensity(f)%ptr, sele)
          oldtdensity_ele_bdy(f,:)=face_val(oldtdensity(f)%ptr, sele)

        end do

        ! loop over the nodes on this surface element
        surface_nodal_loop_i: do iloc = 1, tfield(1)%ptr%mesh%faces%shape%loc

          ! loop over the faces in this surface element
          surface_face_loop: do face = 1, cvfaces%sfaces

            ! is this face a neighbour of iloc?
            if(cvfaces%sneiloc(iloc,face)/=0) then

              ! loop over the gauss pts on this face
              surface_quadrature_loop: do gi = 1, cvfaces%shape%ngi

                ! global gauss point index
                ggi = (face-1)*cvfaces%shape%ngi + gi

                ! u.n
                udotn_bdy=dot_product(u_bdy_f(:,ggi), normal_bdy(:,ggi))

                if(udotn_bdy>0) then
                  income=0.0 ! flow leaving the domain
                else
                  income=1.0 ! flow entering the domain
                end if

                ! as we're on the boundary it's not possible to use high order methods so just
                ! default to the pivotted solution method (first order upwinding)
                ! if the flow is incoming then use the bc ghost values
                ! if the flow is outgoing then use the surface nodes value

                surface_field_loop: do f = 1, nfields
                  
                  if((tfield_bc_type(f,sele)==3)) then
                    ! zero_flux
                    udotn = 0.0
                  else
                    udotn=udotn_bdy
                  end if

                  ! for tfield
                  tfield_face_val(f) = income*ghost_tfield_ele_bdy(f,iloc) + (1.-income)*tfield_ele_bdy(f,iloc)
                  oldtfield_face_val(f) = income*ghost_oldtfield_ele_bdy(f,iloc) + (1.-income)*oldtfield_ele_bdy(f,iloc)

                  ! for tdensity
                  tdensity_face_val = income*ghost_tdensity_ele_bdy(f,iloc) + (1.-income)*tdensity_ele_bdy(f,iloc)
                  oldtdensity_face_val = income*ghost_oldtdensity_ele_bdy(f,iloc) + (1.-income)*oldtdensity_ele_bdy(f,iloc)

                  tdensity_theta_val = tdensity_options(f)%theta*tdensity_face_val + (1.-tdensity_options(f)%theta)*oldtdensity_face_val

                  ! assemble matrix
                  if(getmat(f)) then
                    call addto_diag(A_m(f), nodes_bdy(iloc), &
                                      ptheta(f)*detwei_bdy(ggi)*udotn*(1.-income)*tdensity_theta_val &  ! if iloc is the donor we can do this implicitly
                                    - ptheta(f)*(1.-beta(f))*detwei_bdy(ggi)*udotn_bdy*tdensity_theta_val)
                  end if

                  ! assemble rhs
                  call addto(rhs(f), nodes_bdy(iloc), &
                              -ptheta(f)*udotn*detwei_bdy(ggi)*income*tdensity_theta_val*ghost_tfield_ele_bdy(f,iloc) & ! but we can't if it's the downwind
                              -(1.-ptheta(f))*udotn*detwei_bdy(ggi)*tdensity_theta_val*oldtfield_face_val(f) &
                              +(1.-ptheta(f))*(1.-beta(f))*udotn_bdy*detwei_bdy(ggi)*tdensity_theta_val*oldtfield_ele_bdy(f,iloc))

                end do surface_field_loop

              end do surface_quadrature_loop

            end if ! sneiloc

          end do surface_face_loop

        end do surface_nodal_loop_i

      end do surface_element_loop

      deallocate(x_ele_bdy, detwei_bdy, normal_bdy, u_bdy_f)
      deallocate(nodes_bdy)
      deallocate(tdensity_ele_bdy, oldtdensity_ele_bdy, tfield_ele_bdy, oldtfield_ele_bdy)
      deallocate(ghost_tdensity_ele_bdy, ghost_oldtdensity_ele_bdy, &
                  ghost_tfield_ele_bdy, ghost_oldtfield_ele_bdy)

      deallocate(tfield_bc_type, tdensity_bc_type)
      do f = 1, nfields
        call deallocate(tfield_bc(f))
        call deallocate(tdensity_bc(f))

        call deallocate(tdensity_upwind(f))
        call deallocate(oldtdensity_upwind(f))

        call deallocate(tfield_upwind(f))
        call deallocate(oldtfield_upwind(f))
      end do

      deallocate(x_ele, x_f, detwei, normal, normgi, u_f)
      deallocate(cfl_ele, tfield_ele, oldtfield_ele, tdensity_ele, oldtdensity_ele)
      deallocate(notvisited)

    end subroutine assemble_coupled_advection_m_cv

    subroutine assemble_bassirebay_diffusion_m_cv(D_m, diff_rhs, &
                                       div_m, grad_rhs, &
                                       diffusivity, q_lumpedmass)

      type(csr_matrix), intent(inout) :: D_m
      type(scalar_field), intent(inout) :: diff_rhs

      type(block_csr_matrix), intent(inout) :: div_m
      type(vector_field), intent(inout) :: grad_rhs

      type(tensor_field), intent(in) :: diffusivity
      type(scalar_field), intent(in) :: q_lumpedmass

      logical :: isotropic

      ewrite(1,*) 'in assemble_bassirebay_diffusion_m_cv'

      ! an optimisation that reduces the number of matrix multiplies if we're isotropic
      isotropic=isotropic_field(diffusivity)

      call mult_div_tensorinvscalar_div_T(D_m, div_m, diffusivity, q_lumpedmass, div_m, &
                                          isotropic)

      call mult_div_tensorinvscalar_vector(diff_rhs, div_m, diffusivity, q_lumpedmass, grad_rhs, &
                                           isotropic)

    end subroutine assemble_bassirebay_diffusion_m_cv

    subroutine assemble_auxiliary_m_cv(grad_m_t, grad_rhs, &
                                       tfield, &
                                       cvfaces, cvshape, cvbdyshape, &
                                       x)
      !!< Assembles the matrix and rhs for the auxilliary gradient equation
      !!< for control volumes

      ! inputs/outputs:
      ! the gradient matrix transposed
      type(block_csr_matrix), intent(inout) :: grad_m_t
      ! the rhs of the auxilliary gradient equation
      type(vector_field), intent(inout) :: grad_rhs

      ! the field being solved for
      type(scalar_field), intent(inout) :: tfield

      ! information about cv faces
      type(cv_faces_type), intent(in) :: cvfaces
      ! shape functions for region and surface
      type(element_type), intent(in) :: cvshape, cvbdyshape
      ! the coordinates
      type(vector_field), intent(in) :: x

      ! local memory:
      ! allocatable memory for coordinates, velocity, normals, determinants, nodes
      ! and the cfl number at the gauss pts and nodes
      real, dimension(:,:), allocatable :: x_ele, x_ele_bdy
      real, dimension(:,:), allocatable :: x_f
      real, dimension(:,:), allocatable :: normal, normal_bdy
      real, dimension(:), allocatable :: detwei, detwei_bdy
      real, dimension(:), allocatable :: normgi
      integer, dimension(:), pointer :: nodes
      integer, dimension(:), allocatable :: nodes_bdy
      real, dimension(:), allocatable :: ghost_val

      ! logical array indicating if a face has already been visited by the opposing node
      logical, dimension(:), allocatable :: notvisited

      ! loop integers
      integer :: ele, sele, iloc, oloc, face, gi, ggi, dim

      ! the type of the bc if integrating over domain boundaries
      integer, dimension(:), allocatable :: tfield_bc_type
      ! fields for the bcs over the entire surface mesh
      type(scalar_field) :: tfield_bc

      character(len=OPTION_PATH_LEN) :: l_option_path

      ewrite(1,*) 'in assemble_auxiliary_m_cv'

      l_option_path = trim(tfield%option_path)

      ! allocate memory for assembly
      allocate(x_ele(x%dim,ele_loc(x,1)), &
               x_f(x%dim, cvshape%ngi), &
               detwei(cvshape%ngi), &
               normal(x%dim, cvshape%ngi), &
               normgi(x%dim))
      allocate(notvisited(cvshape%ngi))

      ! Clear memory of arrays being designed
      call zero(grad_m_t)
      call zero(grad_rhs)

      ! loop over elements
      element_loop: do ele=1, element_count(tfield)
        x_ele=ele_val(x, ele)
        x_f=ele_val_at_quad(x, ele, cvshape)
        nodes=>ele_nodes(tfield, ele)

        ! find determinant and unorientated normal
        call transform_cvsurf_to_physical(x_ele, cvshape, &
                                          detwei, normal, cvfaces)

        notvisited=.true.

        ! loop over nodes within this element
        nodal_loop_i: do iloc = 1, tfield%mesh%shape%loc

          ! loop over cv faces internal to this element
          face_loop: do face = 1, cvfaces%faces

            ! is this a face neighbouring iloc?
            if(cvfaces%neiloc(iloc, face) /= 0) then
              oloc = cvfaces%neiloc(iloc, face)

              ! loop over gauss points on face
              quadrature_loop: do gi = 1, cvfaces%shape%ngi

                ! global gauss pt index
                ggi = (face-1)*cvfaces%shape%ngi + gi

                ! have we been here before?
                if(notvisited(ggi)) then
                  notvisited(ggi)=.false.

                  ! correct the orientation of the normal so it points away from iloc
                  normgi=orientate_cvsurf_normgi(x_ele(:,iloc),x_f(:,ggi),normal(:,ggi))

                  ! assemble the matrix
                  dimension_loop1: do dim = 1, mesh_dim(tfield)


                    call addto_diag(grad_m_t, 1, dim, nodes(iloc), &
                              0.5*detwei(ggi)*normgi(dim))
                    call addto(grad_m_t, 1, dim, nodes(oloc), nodes(iloc), &
                              0.5*detwei(ggi)*normgi(dim)) ! remember this is transposed

                    ! notvisited
                    call addto_diag(grad_m_t, 1, dim, nodes(oloc), &
                              -0.5*detwei(ggi)*normgi(dim))
                    call addto(grad_m_t, 1, dim, nodes(iloc), nodes(oloc), &
                              -0.5*detwei(ggi)*normgi(dim)) ! remember this is transposed
                    ! end of notvisited

                  end do dimension_loop1

                end if ! notvisited
              end do quadrature_loop
            end if ! neiloc
          end do face_loop
        end do nodal_loop_i
      end do element_loop

      ! allocate memory for assembly
      allocate(x_ele_bdy(x%dim, face_loc(x,1)), &
              detwei_bdy(cvbdyshape%ngi), &
              normal_bdy(x%dim, cvbdyshape%ngi))
      allocate(tfield_bc_type(surface_element_count(tfield)), &
              nodes_bdy(face_loc(tfield, 1)), &
              ghost_val(face_loc(tfield, 1)))

      ! get the fields over the surface containing the bcs
      call get_entire_boundary_condition(tfield, (/ &
        "weakdirichlet", &
        "neumann      "/), tfield_bc, tfield_bc_type)

      ! loop over the surface elements
      do sele = 1, surface_element_count(tfield)

        ele = face_ele(x, sele)
        x_ele = ele_val(x, ele)
        x_ele_bdy = face_val(x, sele)
        nodes_bdy = face_global_nodes(tfield, sele)
        ghost_val = ele_val(tfield_bc, sele)

        ! calculate the determinant and orientated normal
        call transform_cvsurf_facet_to_physical(x_ele, x_ele_bdy, &
                              cvbdyshape, normal_bdy, detwei_bdy)

        ! loop over the nodes on this surface element
        do iloc = 1, tfield%mesh%faces%shape%loc

          ! loop over the faces in this surface element
          do face = 1, cvfaces%sfaces

            ! is this face a neighbour of iloc?
            if(cvfaces%sneiloc(iloc,face)/=0) then

              ! loop over the gauss pts on this face
              do gi = 1, cvfaces%shape%ngi

                ! global gauss point index
                ggi = (face-1)*cvfaces%shape%ngi + gi

                if(tfield_bc_type(sele)==1) then
                  ! assemble grad_rhs
                  call addto(grad_rhs, nodes_bdy(iloc), &
                              detwei_bdy(ggi)*normal_bdy(:,ggi)*ghost_val(iloc))



                else

                  dimension_loop2: do dim = 1, mesh_dim(tfield)
                    ! assemble matrix
                    call addto_diag(grad_m_t, 1, dim, nodes_bdy(iloc), &
                                      detwei_bdy(ggi)*normal_bdy(dim,ggi))
                  end do dimension_loop2

                end if

              end do ! gi

            end if ! sneiloc

          end do ! face

        end do ! iloc

      end do ! sele

      deallocate(x_ele_bdy, detwei_bdy, normal_bdy)
      deallocate(nodes_bdy)

      deallocate(tfield_bc_type)
      call deallocate(tfield_bc)

      deallocate(x_ele, x_f, detwei, normal, normgi)
      deallocate(notvisited)

    end subroutine assemble_auxiliary_m_cv
    ! end of assembly subroutines
    !************************************************************************
    subroutine calculate_auxiliary_gradient(state, aux_grad)

      type(state_type), intent(inout) :: state
      type(vector_field), intent(inout) :: aux_grad

      type(block_csr_matrix) :: grad_m_t
      type(vector_field) :: grad_rhs

      character(len=FIELD_NAME_LEN) :: field_name
      type(scalar_field), pointer :: sfield
      type(csr_sparsity), pointer :: mesh_sparsity
      type(vector_field), pointer :: x

      integer :: quaddegree, dim
      type(cv_faces_type) :: cvfaces
      type(element_type) :: cvshape, cvbdyshape

      type(scalar_field) :: lumpedmass

      type(scalar_field) :: tempfield, temprhs

      call get_option(trim(aux_grad%option_path)//"/diagnostic/gradient_of_field", &
                      field_name)

      sfield => extract_scalar_field(state, trim(field_name))
      mesh_sparsity=>get_csr_sparsity_firstorder(state, sfield%mesh, sfield%mesh)
      x => extract_vector_field(state, "Coordinate")

      call allocate(grad_m_t, sparsity=mesh_sparsity, blocks=(/1,mesh_dim(sfield)/), &
                    name=trim(field_name)//"AuxiliaryGradientMatrixTransposed")
      call zero(grad_m_t)

      call allocate(grad_rhs, mesh_dim(sfield), sfield%mesh, &
                    name=trim(field_name)//"AuxiliaryGradientRHS")
      call zero(grad_rhs)

      call allocate(lumpedmass, sfield%mesh, name=trim(field_name)//"LumpedMass")
      call zero(lumpedmass)

      call get_option("/geometry/quadrature/controlvolume_surface_degree", &
                     quaddegree, default=1)

      cvfaces=find_cv_faces(vertices=ele_vertices(sfield, 1), &
                            dimension=sfield%mesh%shape%dim, &
                            polydegree=sfield%mesh%shape%degree, &
                            quaddegree=quaddegree)

      cvshape=make_cv_element_shape(cvfaces, x%mesh%shape%degree)
      cvbdyshape=make_cvbdy_element_shape(cvfaces, x%mesh%faces%shape%degree)

      call compute_lumped_mass(x, lumpedmass)

      call assemble_auxiliary_m_cv(grad_m_t, grad_rhs, &
                                   sfield, &
                                   cvfaces, cvshape, cvbdyshape, &
                                   x)

      ewrite_minmax(sfield%val)

      do dim = 1, aux_grad%dim
        tempfield = extract_scalar_field(aux_grad, dim)
        temprhs = extract_scalar_field(grad_rhs, dim)
        call mult_T(tempfield, block(grad_m_t,1,dim), sfield)
        call addto(tempfield, temprhs)
        tempfield%val = tempfield%val/lumpedmass%val
      end do

      call deallocate(grad_m_t)
      call deallocate(grad_rhs)
      call deallocate(lumpedmass)
      call deallocate(cvshape)
      call deallocate(cvbdyshape)
      call deallocate(cvfaces)

    end subroutine calculate_auxiliary_gradient
    !************************************************************************
    ! subroutines dealing with the writing of the advection_convergence files
    subroutine initialise_advection_convergence(state)

      type(state_type), dimension(:), intent(in) :: state

      integer :: nfiles

      logical, save :: initialised=.false.

      integer :: column, i, j, fileno
      character(len=254) :: buffer
      character(len=FIELD_NAME_LEN) :: material_phase_name, field_name

      if(initialised) return
      initialised=.true.

      nfiles = option_count("/material_phase/scalar_field/prognostic/output/convergence_file")

      allocate(conv_unit(nfiles))
      allocate(sfield_list(nfiles))

      if(nfiles==0) return

      fileno = 0
      do i = 1, size(state)
        material_phase_name=trim(state(i)%name)

        do j = 1, size(state(i)%scalar_fields)

          if(have_option(trim(state(i)%scalar_fields(j)%ptr%option_path)//&
              "/prognostic/output/convergence_file")) then
            field_name=trim(state(i)%scalar_fields(j)%ptr%name)

            fileno=fileno+1

            if(fileno>nfiles) then
              ewrite(-1,*) 'fileno = ', fileno, 'nfiles = ', nfiles
              FLAbort("More fields think they want a convergence file than expected.")
            end if

            sfield_list(fileno) = trim(material_phase_name)//&
                      "__"//trim(field_name)

            ! open and write a file (if its the first processor)
            if(getprocno() == 1) then
              conv_unit(fileno) = free_unit()
              open(unit=conv_unit(fileno), file=trim(sfield_list(fileno))//".convergence", &
                      action="write")


              write(conv_unit(fileno), '(a)') "<header>"

              column=0
              ! Initial columns are elapsed time, dt, global iteration and advective iteration
              column=column+1
              buffer=field_tag(name="ElapsedTime", column=column, statistic="value")
              write(conv_unit(fileno), '(a)') trim(buffer)
              column=column+1
              buffer=field_tag(name="dt", column=column, statistic="value")
              write(conv_unit(fileno), '(a)') trim(buffer)
              column=column+1
              buffer=field_tag(name="Iteration", column=column, statistic="value")
              write(conv_unit(fileno), '(a)') trim(buffer)
              column=column+1
              buffer=field_tag(name="Subcycle", column=column, statistic="value")
              write(conv_unit(fileno), '(a)') trim(buffer)
              column=column+1
              buffer=field_tag(name="AdvectionIteration", column=column, statistic="value")
              write(conv_unit(fileno), '(a)') trim(buffer)

              column=column+1
              buffer=field_tag(name=trim(field_name), column=column, statistic="error", material_phase_name=trim(material_phase_name))
              write(conv_unit(fileno), '(a)') trim(buffer)

              write(conv_unit(fileno), '(a)') "</header>"

            end if

          end if
        end do
      end do

      if(fileno/=nfiles) then
        ewrite(-1,*) 'fileno = ', fileno, 'nfiles = ', nfiles
        FLAbort("Fewer fields thought they wanted a convergence file than expected.")
      end if

    end subroutine initialise_advection_convergence

    subroutine test_and_write_advection_convergence(field, nlfield, coordinates, filename, &
                                                    time, dt, it, subcyc, adv_it, &
                                                    error)

       type(scalar_field), intent(inout) :: field, nlfield
       type(vector_field), intent(in) :: coordinates
       character(len=*), intent(in) :: filename
       real, intent(in) :: time, dt
       integer, intent(in) :: it, subcyc, adv_it

       real, intent(out) :: error

       logical :: write_convergence_file
       character(len=10) :: format, iformat
       integer :: fileno
       
       integer :: convergence_norm
       
       convergence_norm = convergence_norm_integer(trim(field%option_path)//&
                          "/prognostic/temporal_discretisation/control_volumes/number_advection_iterations/tolerance")

       error = 0.0
       call field_con_stats(field, nlfield, error, &
                            convergence_norm, coordinates)

       format='(e15.6e3)'
       iformat='(i4)'

       write_convergence_file = .false.
       fileno=find_fileno(filename)
       if(fileno/=0) then
         write_convergence_file = .true.
       end if

       if(write_convergence_file) then
         if(getprocno() == 1) then
           write(conv_unit(fileno), format, advance="no") time
           write(conv_unit(fileno), format, advance="no") dt
           write(conv_unit(fileno), iformat, advance="no") it
           write(conv_unit(fileno), iformat, advance="no") subcyc
           write(conv_unit(fileno), iformat, advance="no") adv_it
           write(conv_unit(fileno), format, advance="no") error
           write(conv_unit(fileno),'(a)') ""  ! end of line
         end if
       end if

    end subroutine test_and_write_advection_convergence

    pure function find_fileno(filename) result(fileno)

      integer :: fileno
      character(len=*), intent(in) :: filename

      integer :: i

      fileno = 0

      do i = 1, size(sfield_list)
        if(trim(filename)==trim(sfield_list(i))) then
          fileno = i
          return
        end if
      end do

    end function find_fileno
    ! end of convergence file subroutines
    !************************************************************************
    !************************************************************************
    ! control volume options checking
    subroutine field_equations_cv_check_options
      integer :: nmat, nfield, m, f, stat
      character(len=OPTION_PATH_LEN) :: mat_name, field_name, diff_scheme
      integer :: equation_type, weakdirichlet_count
      logical :: cv_disc, mmat_cv_disc, diff, conv_file, subcycle, cv_temp_disc, tolerance, explicit
      real :: theta, p_theta

      nmat = option_count("/material_phase")

      do m = 0, nmat-1
         call get_option("/material_phase["//int2str(m)//"]/name", mat_name)
         nfield = option_count("/material_phase["//int2str(m)//"]/scalar_field")
         do f = 0, nfield-1
            call get_option("/material_phase["//int2str(m)//"]/scalar_field["//int2str(f)//&
                            "]/name", field_name)
            cv_disc=have_option("/material_phase["//int2str(m)//"]/scalar_field["//int2str(f)//&
                            "]/prognostic/spatial_discretisation/control_volumes")
            mmat_cv_disc=have_option("/material_phase["//int2str(m)//"]/scalar_field["//int2str(f)//&
                            "]/prognostic/spatial_discretisation/coupled_cv")
            diff=have_option("/material_phase["//int2str(m)//"]/scalar_field["//int2str(f)//&
                            "]/prognostic/tensor_field::Diffusivity")
            call get_option("/material_phase["//int2str(m)//"]/scalar_field["//int2str(f)//&
                            "]/prognostic/spatial_discretisation/control_volumes/diffusion_scheme[0]/name", &
                            diff_scheme, default="None")
            conv_file=have_option("/material_phase["//int2str(m)//"]/scalar_field["//int2str(f)//&
                            "]/prognostic/output/convergence_file")

            cv_temp_disc=have_option("/material_phase["//int2str(m)//"]/scalar_field["//int2str(f)//&
                            "]/prognostic/temporal_discretisation/control_volumes")
            tolerance=have_option("/material_phase["//int2str(m)//"]/scalar_field["//int2str(f)//&
                            "]/prognostic/temporal_discretisation/control_volumes/number_advection_iterations/tolerance")
            subcycle=((have_option("/material_phase["//int2str(m)//"]/scalar_field["//int2str(f)//&
                            "]/prognostic/temporal_discretisation/control_volumes/maximum_courant_number_per_subcycle")).or.&
                      (have_option("/material_phase["//int2str(m)//"]/scalar_field["//int2str(f)//&
                            "]/prognostic/temporal_discretisation/control_volumes/number_advection_subcycles")))
            explicit=have_option("/material_phase["//int2str(m)//"]/scalar_field["//int2str(f)//&
                            "]/prognostic/explicit")

            weakdirichlet_count=option_count("/material_phase["//int2str(m)//"]/scalar_field["//int2str(f)//&
                            "]/prognostic/boundary_conditions/type[0]/apply_weakly")

            equation_type=equation_type_index(trim("/material_phase["//int2str(m)//"]/scalar_field["//int2str(f)//"]"))
            select case(equation_type)
            case(FIELD_EQUATION_CONSERVATIONOFMASS, FIELD_EQUATION_REDUCEDCONSERVATIONOFMASS, FIELD_EQUATION_INTERNALENERGY)
              if(cv_disc.and.diff) then
                ewrite(-1,*) "Options checking field "//&
                              trim(field_name)//" in material_phase "//&
                              trim(mat_name)//"."
                FLExit("Selected equation type not compatible with diffusion")
              end if
            end select

            if(mmat_cv_disc) then
              if(diff) then
                ewrite(-1,*) "Options checking field "//&
                              trim(field_name)//" in material_phase "//&
                              trim(mat_name)//"."
                ewrite(-1,*) "Use control volume discretisation if you want diffusion."
                FLExit("Multiple coupled control volume discretisation not yet compatible with Diffusivity")
              end if
              
              if(.not.have_option("/material_phase["//int2str(m)//"]/scalar_field["//int2str(f)//&
                            "]/prognostic/priority")) then
                FLExit("Coupled control volume discretisation requires a priority option.")
              end if
              
              if(explicit) then
                
                call get_option("/material_phase["//int2str(m)//"]/scalar_field["//int2str(f)//&
                            "]/prognostic/temporal_discretisation/theta", theta)
                if (theta/=0.0) then
                  ewrite(-1,*) "Options checking field "//&
                                trim(field_name)//" in material_phase "//&
                                trim(mat_name)//"."
                  FLExit("Explicit coupled control volume discretisations must use temporal_discretisation/theta = 0.0")
                end if
                
                call get_option("/material_phase["//int2str(m)//"]/scalar_field["//int2str(f)//&
                            "]/prognostic/temporal_discretisation/control_volumes/pivot_theta", p_theta, default=1.0)
                if (p_theta/=0.0) then
                  ewrite(-1,*) "Options checking field "//&
                  trim(field_name)//" in material_phase "//&
                  trim(mat_name)//"."
                  FLExit("Explicit coupled control volume discretisations must use temporal_discretisation/control_volumes/pivot_theta = 0.0")
                end if
                
              end if
            elseif(cv_disc) then
              if(diff) then
                select case(diff_scheme)
                case("ElementGradient")
                  if(weakdirichlet_count>0) then
                    ewrite(-1,*) "Options checking field "//&
                                  trim(field_name)//" in material_phase "//&
                                  trim(mat_name)//"."
                    ewrite(-1,*) "ElementGradient diffusion scheme not compatible with weak dirichlet boundary conditions!"
                    ewrite(-1,*) "Use strong dirichlet boundary conditions or switch the diffusion scheme to BassiRebay."
                    ewrite(-1,*) "Sorry and Good Luck!"
                    FLExit("ElementGradient diffusion scheme not compatible with weak dirichlet boundary conditions")
                  end if
                end select
              end if
              if(explicit) then
                
                call get_option("/material_phase["//int2str(m)//"]/scalar_field["//int2str(f)//&
                            "]/prognostic/temporal_discretisation/theta", theta)
                if (theta/=0.0) then
                  ewrite(-1,*) "Options checking field "//&
                                trim(field_name)//" in material_phase "//&
                                trim(mat_name)//"."
                  FLExit("Explicit control volume discretisations must use temporal_discretisation/theta = 0.0")
                end if
                
                call get_option("/material_phase["//int2str(m)//"]/scalar_field["//int2str(f)//&
                            "]/prognostic/temporal_discretisation/control_volumes/pivot_theta", p_theta, default=1.0)
                if (p_theta/=0.0) then
                  ewrite(-1,*) "Options checking field "//&
                  trim(field_name)//" in material_phase "//&
                  trim(mat_name)//"."
                  FLExit("Explicit control volume discretisations must use temporal_discretisation/control_volumes/pivot_theta = 0.0")
                end if
                
              end if
            else
              if(conv_file) then
                ewrite(-1,*) "Options checking field "//&
                              trim(field_name)//" in material_phase "//&
                              trim(mat_name)//"."
                FLExit("Only pure control volume and coupled_cv discretisations can output a convergence file")
              end if
              if(cv_temp_disc) then
                ewrite(-1,*) "Options checking field "//&
                              trim(field_name)//" in material_phase "//&
                              trim(mat_name)//"."
                FLExit("Only control volume or coupled_cv discretisations can use control_volume temporal discretisations")
              end if
              if(explicit) then
                ewrite(-1,*) "Options checking field "//&
                              trim(field_name)//" in material_phase "//&
                              trim(mat_name)//"."
                FLExit("Only pure control volume or coupled_cv discretisations can solve explicitly")
              end if
            end if

         end do
       end do

    end subroutine field_equations_cv_check_options
    ! end of control volume options checking
    !************************************************************************

end module field_equations_CV
