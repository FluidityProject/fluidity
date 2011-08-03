
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


module multiphase_field_advection

  use multiphase_1D_engine
  use printout

  use fldebug
  use state_module
  use fields
  use spud
  use global_parameters, only: option_path_len
  use advection_diffusion_cg
  use advection_diffusion_DG
  use advection_diffusion_FV
  use field_equations_cv, only: solve_field_eqn_cv

  implicit none

  private 

  public  :: solve_multiphase_field_advection

contains

  subroutine solve_multiphase_field_advection( state, problem, nphase, ncomp, totele, ndim, &
       u_nloc, xu_nloc, cv_nloc, x_nloc, mat_nloc, &
       cv_snloc, u_snloc, stotel, &
       domain_length, &
                                ! Element types
       u_ele_type, p_ele_type, cv_ele_type, &
       cv_sele_type, u_sele_type, &
                                ! Total time loop and initialisation parameters
       ntime, ntime_dump, nits, &
       nits_flux_lim_volfra, nits_flux_lim_comp, &
       dt,  &
                                ! Discretisation parameters
       t_beta, t_theta, t_disopt, t_dg_vel_int_opt, lump_eqns, &
       use_theta_flux, get_theta_flux, &
       opt_vel_upwind_coefs, nopt_vel_upwind_coefs, &
       noit_dim, &
       sat_error_relax2_noit, t_error_relax2_noit, gl_error_relax2_noit, &
       u_error_relax2_noit, p_error_relax2_noit, mass_error_relax2_noit, &
       in_ele_upwind, dg_ele_upwind, &
                                ! Total nodes for different meshes
       cv_nonods, u_nonods, mat_nonods, x_nonods, &
       u_ndgln, xu_ndgln, cv_ndgln, x_ndgln, &
       mat_ndgln, u_sndgln, cv_sndgln, &
                                ! Boundary conditions and surface elements
       wic_d_bc, wic_u_bc, wic_t_bc, &
       suf_d_bc, suf_t_bc, suf_u_bc, suf_v_bc, suf_w_bc, &
       suf_t_bc_rob1, suf_t_bc_rob2, &
                                ! Positions and grid velocities
       x, y, z, nu, nv, nw, ug, vg, wg, &
                                ! Absorption and source terms and coefficients
       t_absorb, t_source, &
                                ! Diffusion parameters
       tdiffusion, &
                                ! Scalar fields
       t, p, cv_one, volfra_pore, deriv, &
       told, nuold, nvold, nwold, &
                                ! Matrices sparsity
       mx_ncolacv, ncolacv, finacv, colacv, midacv, & ! CV multi-phase eqns (e.g. vol frac, temp)
       mx_nct, ncolct, findct, colct, & ! CT sparsity - global cty eqn
       ncolm, findm, colm, midm, & ! CV-FEM matrix
       mxnele, ncolele, finele, colele, & ! Element connectivity 
       option_path)
       
    implicit none
    
    !! New variable declaration

    type(state_type), dimension(:), intent(inout) :: state

    !! Old variable declaration

    integer, intent( in ) :: problem, nphase, ncomp, totele, ndim, u_nloc, xu_nloc, cv_nloc, x_nloc, &
         mat_nloc, cv_snloc, u_snloc, stotel
    real, intent( in ) :: domain_length
    integer, intent( in ) :: u_ele_type, p_ele_type, cv_ele_type, cv_sele_type, u_sele_type, &
         ntime, ntime_dump, nits, nits_flux_lim_volfra, nits_flux_lim_comp
    real :: dt
    ! The following need to be changed later in the other subrts as it should be controlled by the 
    ! input files
    real, intent( inout ) :: t_beta, t_theta
    integer, intent( in ) :: t_disopt, t_dg_vel_int_opt
    logical, intent( in ) :: lump_eqns, use_theta_flux, get_theta_flux
    integer, intent( in ) :: nopt_vel_upwind_coefs, cv_nonods, u_nonods, mat_nonods, x_nonods
    real, dimension( nopt_vel_upwind_coefs ), intent( inout ) :: opt_vel_upwind_coefs
    integer, intent( in ) :: noit_dim
    real, dimension( noit_dim ), intent( in ) :: sat_error_relax2_noit, t_error_relax2_noit, gl_error_relax2_noit, &
         u_error_relax2_noit, p_error_relax2_noit, mass_error_relax2_noit
    integer, intent( in ) :: in_ele_upwind, dg_ele_upwind

    integer, dimension( totele * u_nloc ), intent( in )  :: u_ndgln
    integer, dimension( totele * xu_nloc ), intent( in )  :: xu_ndgln
    integer, dimension( totele * cv_nloc ), intent( in )  :: cv_ndgln
    integer, dimension( totele * cv_nloc ), intent( in )  :: x_ndgln
    integer, dimension( totele * mat_nloc ), intent( in )  :: mat_ndgln
    integer, dimension( totele * u_snloc ), intent( in )  :: u_sndgln
    integer, dimension( totele * cv_snloc ), intent( in )  :: cv_sndgln
    integer, dimension( stotel * nphase ), intent( in ) :: wic_d_bc, wic_u_bc, wic_t_bc
    real, dimension( stotel * cv_snloc * nphase ), intent( in ) :: suf_d_bc, suf_t_bc
    real, dimension( stotel * u_snloc * nphase ), intent( in ) :: suf_u_bc, suf_v_bc, suf_w_bc
    real, dimension( stotel * cv_snloc * nphase ), intent( in ) :: suf_t_bc_rob1, suf_t_bc_rob2
    ! Variables bellow may change on this subroutine, however this should be changed later
    real, dimension( x_nonods ), intent( inout ) :: x, y, z
    real, dimension( u_nonods * nphase ), intent( inout ) :: nu, nv, nw
    real, dimension( cv_nonods, nphase, nphase ), intent( inout ) :: t_absorb
    real, dimension( cv_nonods * nphase ), intent( inout ) :: t_source
    real, dimension( mat_nonods, ndim, ndim, nphase ), intent( inout ) :: tdiffusion
    ! Variables initialised before but will change during computational time
    real, dimension( cv_nonods * nphase ), intent( inout ) :: t, told, cv_one
    real, dimension( cv_nonods ), intent( inout ) :: p
    real, dimension( totele ), intent( inout ) :: volfra_pore
    real, dimension( cv_nonods * nphase ), intent( inout ) :: deriv
    real, dimension( u_nonods * nphase ), intent( inout ) :: nuold, nvold, nwold, ug, vg, wg
    integer, intent ( in ) :: mx_ncolacv, ncolacv
    integer, dimension( cv_nonods * nphase + 1 ), intent (in ) :: finacv
    integer, dimension( mx_ncolacv ), intent (in ) :: colacv
    integer, dimension( cv_nonods * nphase ), intent (in ) :: midacv
    integer, intent ( in ) :: mx_nct, ncolct
    integer, dimension( cv_nonods + 1 ), intent (in ) :: findct
    integer, dimension( mx_nct ), intent (in ) :: colct
    integer, intent ( in ) :: ncolm
    integer, dimension( cv_nonods + 1 ), intent (in ) :: findm
    integer, dimension( ncolm ), intent (in ) :: colm
    integer, dimension( cv_nonods ), intent (in ) :: midm
    integer, intent ( in ) :: mxnele, ncolele
    integer, dimension( totele + 1 ), intent (in ) :: finele
    integer, dimension( mxnele ), intent (in ) :: colele
    character(len=*), intent(in), optional :: option_path
    
    ! Local variables
    real :: acctim, xacc
    integer :: its, itime, ele, iloc, x_nod, cv_nod
    character( len = 100 ) :: field
    integer :: unit_T
    real, dimension( : ), allocatable :: T_FEMT, MEAN_PORE_CV

    REAL, DIMENSION( :, :, :, : ), allocatable :: THETA_FLUX, ONE_M_THETA_FLUX
    INTEGER :: IGOT_T2, IGOT_THETA_FLUX, SCVNGI_THETA
    
    logical :: do_fluidity_advection
    real :: finish_time

    IGOT_T2 = 0
    IGOT_THETA_FLUX = 0 
    SCVNGI_THETA = 0

    ALLOCATE( THETA_FLUX( TOTELE * IGOT_THETA_FLUX, CV_NLOC, SCVNGI_THETA, NPHASE ))
    ALLOCATE( ONE_M_THETA_FLUX( TOTELE * IGOT_THETA_FLUX, CV_NLOC, SCVNGI_THETA, NPHASE ))
    ALLOCATE( T_FEMT(  cv_nonods * nphase ))
    ALLOCATE( MEAN_PORE_CV(  cv_nonods ))


    call get_option("/timestepping/current_time", acctim)
    call get_option("/timestepping/timestep", dt)
    call get_option("/timestepping/finish_time", finish_time)
    
    ewrite(3,*) 'In solve_multiphase_field_advection'
    
    itime = 0
    
    Loop_Time: do 
       
       itime = itime + 1
       
       acctim = acctim + dt
       
       ewrite(3,*) '********************************'
       ewrite(3,*) 'NEW TIMESTEP'
       ewrite(3,*) 'acctim: ', acctim
       ewrite(3,*) 'finish_time: ', finish_time
       ewrite(3,*) '********************************'
       
       if (acctim > finish_time) exit Loop_Time 
       
       told = t

       Loop_Non_Linear_Iterations: do its = 1, nits
          write( 375, * ) 'itime, its: ', itime, its
          
          !! Put alternative call to solve_advection_diffusion_? or
          !! solve_field_eqn_c? as in Fluids.F90
          !! Will need to pull out and put back fields at each timestep.
          !!
          !! Hopefully this approach will be portable to some extent to the
          !! call to intenerge_assem_solve in Multiphase_Mom_Press_Volf
          do_fluidity_advection = .true.
          
          if (do_fluidity_advection) then
             call use_fluidity_advection_routines(state, t, cv_ndgln, its, dt)
          
          else
             call INTENERGE_ASSEM_SOLVE(  &
               NCOLACV, FINACV, COLACV, MIDACV, & ! CV sparsity pattern matrix
               NCOLCT, FINDCT, COLCT, &
               CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
               U_ELE_TYPE, CV_ELE_TYPE, CV_SELE_TYPE,  &
               NPHASE,  &
               CV_NLOC, U_NLOC, X_NLOC,  &
               CV_NDGLN, X_NDGLN, U_NDGLN, &
               CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
               X, Y, Z, &
               NU, NV, NW, NUOLD, NVOLD, NWOLD, UG, VG, WG, &
               T, TOLD, &
               CV_ONE, CV_ONE, &
               MAT_NLOC, MAT_NDGLN, MAT_NONODS, TDIFFUSION, &
               T_DISOPT, T_DG_VEL_INT_OPT, DT, T_THETA, T_BETA, &
               SUF_T_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, &
               SUF_T_BC_ROB1, SUF_T_BC_ROB2,  &
               WIC_T_BC, WIC_D_BC, WIC_U_BC, &
               DERIV, P, &
               T_SOURCE, T_ABSORB, VOLFRA_PORE, &
               NDIM,  &
               NCOLM, FINDM, COLM, MIDM, &
               XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, LUMP_EQNS, &
               OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, T_FEMT, CV_ONE, &
               IGOT_T2,T,TOLD, IGOT_THETA_FLUX, SCVNGI_THETA, GET_THETA_FLUX, USE_THETA_FLUX, &
               T, T, T, &
               SUF_T_BC, SUF_T_BC_ROB1, SUF_T_BC_ROB2, WIC_T_BC, IN_ELE_UPWIND, DG_ELE_UPWIND, &
               NOIT_DIM, &
               T_ERROR_RELAX2_NOIT, MASS_ERROR_RELAX2_NOIT, NITS_FLUX_LIM_COMP, &
               MEAN_PORE_CV, &
               option_path = trim(option_path) )
          endif

       end do Loop_Non_Linear_Iterations

       Select Case( problem )
       Case( -1 );
          field = 'CV_adv_field_Tdg'
       Case( 0 );
          field = 'CV_adv_field_Tstd'
       Case( -2 );
          field = 'CV_adv_field_Tcty'            
       end Select

       unit_T = 1
       call generate_name_dump( itime, unit_T, field, 1 )
       xacc = 0.
       do ele = 1, totele
          do iloc = 1, cv_nloc
             cv_nod = cv_ndgln(( ele - 1 ) * cv_nloc + iloc )
             x_nod = x_ndgln(( ele - 1 ) * cv_nloc + iloc )
             if( problem == -2 ) then
                write( itime, * ) xacc, t( cv_nod )
                xacc = unit_T + ( domain_length / real( totele ) ) / real( cv_nloc - 1 )
             else            
                write( unit_T, * ) x( x_nod ), t( cv_nod )
             end if
          end do
          xacc = xacc - ( domain_length / real( totele ) ) / real( cv_nloc - 1 )
       end do

    end do Loop_Time
    
    call set_option("/timestepping/current_time", ACCTIM)
    call set_option("/timestepping/timestep", dt)  
      
    ewrite(3,*) 'Leaving solve_multiphase_field_advection'

  end subroutine solve_multiphase_field_advection
  
  subroutine use_fluidity_advection_routines(state, t, cv_ndgln, its, dt)
  
    implicit none
  
    type(state_type), dimension(:), intent(inout) :: state
    type(scalar_field), pointer :: temperature
    integer, intent(in) :: its
    real, intent(in) :: dt
    integer, dimension(:), intent(in) :: cv_ndgln
    real, dimension(:), intent(inout) :: t
    
    !local variables
    integer :: i, j, stat
    integer, dimension(:), pointer :: element_nodes
    
    character(len=OPTION_PATH_LEN) :: option_path, field_name
    
    !! Pull temperature out of state
    temperature => extract_scalar_field(state, "Temperature", stat)
    if (stat/=0) FLAbort("Advected field not called Temperature!")
    
    option_path = "/material_phase[0]/scalar_field::Temperature"
    field_name = "Temperature"
    
    !! Update temperature from prototype field
    temp_ele_loop: do i = 1,element_count(temperature)
            
      element_nodes => ele_nodes(temperature,i)
            
      temp_node_loop: do j = 1,size(element_nodes)
             
        call set(temperature, &
                 element_nodes(j), &
                 t((cv_ndgln((i-1)*size(element_nodes)+j)))) ! + (p-1)*node_count(temperature)))
            
      end do temp_node_loop
            
    end do temp_ele_loop

    !! Call appropriate advection routine
    if(have_option(trim(option_path)//&
                     & "/prognostic/spatial_discretisation/discontinuous_galerkin")) then

      ! Solve the DG form of the equations.
      call solve_advection_diffusion_dg(field_name=trim(field_name), &
                        & state=state(1))

    elseif(have_option(trim(option_path)//&
                     & "/prognostic/spatial_discretisation/finite_volume")) then

      ! Solve the FV form of the equations.
      call solve_advection_diffusion_fv(field_name=trim(field_name), &
                        & state=state(1))

    elseif(have_option(trim(option_path)//&
                     & "/prognostic/spatial_discretisation/control_volumes")) then

      ! Solve the pure control volume form of the equations
      call solve_field_eqn_cv(field_name=trim(field_name), &
                        state=state, &
                        global_it=its)

    elseif(have_option(trim(option_path) // &
                     & "/prognostic/spatial_discretisation/continuous_galerkin")) then

      call solve_field_equation_cg(trim(field_name), state(1), dt)
      
    else
      ewrite(2, *) "Not solving scalar field " // trim(field_name) // " in state " // trim(state(1)%name) //" in an advdif-like subroutine."

    end if ! End of dg/cv/cg choice.
    
    
    !! Update prototype field from temperature
    t_ele_loop: do i = 1,element_count(temperature)
            
      element_nodes => ele_nodes(temperature,i)
            
      t_node_loop: do j = 1,size(element_nodes)
                 
        t((cv_ndgln((i-1)*size(element_nodes)+j))) = & !+ (p-1)*node_count(temperature)) = &
                        temperature%val(element_nodes(j))
            
      end do t_node_loop
            
    end do t_ele_loop
  
  end subroutine use_fluidity_advection_routines

end module multiphase_field_advection
