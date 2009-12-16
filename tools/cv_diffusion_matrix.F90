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

program cv_diffusion_matrix
  use spud
  use populate_state_module
  use global_parameters
  use state_module
  use unittest_tools
  use vtk_interfaces
  use field_equations_cv
  use cv_options
  use sparsity_patterns
  use fldebug
  use mpi_interfaces
  use boundary_conditions
  use cv_shape_functions
  use cv_faces
  use cvtools
  use cv_upwind_values
  use cv_face_values
  use fefields, only: compute_lumped_mass
  use reference_counting
  use sparse_matrices_fields

  implicit none

#ifdef HAVE_MPI
  integer :: ierr
#endif
  type(state_type), dimension(:), pointer :: states => null()

  character(len = 512) :: filename, simname

  type(scalar_field), target :: dummydensity
  type(scalar_field), pointer :: tfield, old_tfield
  type(scalar_field), pointer :: tdensity, old_tdensity

  type(vector_field), pointer :: x
  type(vector_field) :: advu

  type(tensor_field), pointer :: diffusivity

  type(cv_options_type) :: tfield_options, tdensity_options

  ! control volume face information
  type(cv_faces_type) :: cvfaces

  ! control volume shape function for volume and boundary
  type(element_type) :: u_cvshape, u_cvbdyshape
  type(element_type) :: x_cvshape, x_cvbdyshape, &
                        x_cvshape_full, x_cvbdyshape_full
  ! t_cvshape is the element with reduced numbers of derivatives
  ! taken across the control volume faces
  ! t_cvshape_full contains the derivatives with respect to the parent
  ! elements canonical coordinates evaluated at the control volume faces
  ! t_cvbdyshape_full is the same but on the boundary
  type(element_type) :: t_cvshape, t_cvshape_full, t_cvbdyshape_full

  type(scalar_field) :: lumpedmass, q_lumpedmass, rhs, diff_rhs, cfl_no
  type(csr_matrix) :: A_m, D_m, mass, inverse_mass
  type(csr_sparsity) :: mesh_sparsity, mesh_sparsity_1, grad_m_t_sparsity

  integer :: quaddegree

  logical :: diffusion, fdcflno, getmat

#ifdef HAVE_MPI
  call MPI_Init(ierr)
  assert(ierr == MPI_SUCCESS)
#endif
  call python_init()

  call set_global_debug_level(0)

  call read_command_line(filename)

  call load_options(filename)
  call populate_state(states)
  call allocate_and_insert_auxilliary_fields(states)

  call get_option("/simulation_name", simname)

  tfield => extract_scalar_field(states(1), "Tracer")
  old_tfield => extract_scalar_field(states(1), "OldTracer")

  call allocate(dummydensity, tfield%mesh, name="DummyDensity", field_type=FIELD_TYPE_CONSTANT)
  call set(dummydensity, 1.0)
  dummydensity%option_path = " "

  tdensity=>dummydensity
  old_tdensity=>dummydensity

  ! now we can get the options for these fields
  ! handily wrapped in a new type...
  tfield_options=get_cv_options(tfield%option_path, tfield%%mesh%shape%numbering%family)
  tdensity_options=get_cv_options(tdensity%option_path, tdensity%mesh%shape%numbering%family)

  x=>extract_vector_field(states(1), "Coordinate")

  ! find relative velocity
  call allocate(advu, mesh_dim(tfield), tfield%mesh, "DummyRelativeVelocity", field_type=FIELD_TYPE_CONSTANT)
  call zero(advu)
  advu%option_path = ""

  diffusivity=>extract_tensor_field(states(1), "TracerDiffusivity")
  diffusion = .true.

  ! create control volume shape functions
  call get_option("/geometry/quadrature/controlvolume_surface_degree", &
                  quaddegree, default=1)
  cvfaces=find_cv_faces(vertices=ele_vertices(tfield, 1), &
                        dimension=mesh_dim(tfield), &
                        polydegree=tfield%mesh%shape%degree, &
                        quaddegree=quaddegree)

  u_cvshape=make_cv_element_shape(cvfaces, advu%mesh%shape%degree)
  x_cvshape=make_cv_element_shape(cvfaces, x%mesh%shape%degree)
  t_cvshape=make_cv_element_shape(cvfaces, tfield%mesh%shape%degree)

  u_cvbdyshape=make_cvbdy_element_shape(cvfaces, advu%mesh%faces%shape%degree)
  x_cvbdyshape=make_cvbdy_element_shape(cvfaces, x%mesh%faces%shape%degree)

  if(tfield_options%diffusionscheme==CV_DIFFUSION_ELEMENTGRADIENT) then
    x_cvshape_full=make_cv_element_shape(cvfaces, x%mesh%shape%degree, &
                   type=ELEMENT_CONTROLVOLUME_SURFACE_BODYDERIVATIVES)
    t_cvshape_full=make_cv_element_shape(cvfaces, tfield%mesh%shape%degree, &
                   type=ELEMENT_CONTROLVOLUME_SURFACE_BODYDERIVATIVES)

    x_cvbdyshape_full=make_cvbdy_element_shape(cvfaces, x%mesh%faces%shape%degree, &
                      type=ELEMENT_CONTROLVOLUME_SURFACE_BODYDERIVATIVES)
    t_cvbdyshape_full=make_cvbdy_element_shape(cvfaces, tfield%mesh%shape%degree, &
                      type=ELEMENT_CONTROLVOLUME_SURFACE_BODYDERIVATIVES)
  else
    x_cvshape_full=x_cvshape
    t_cvshape_full=t_cvshape
    x_cvbdyshape_full=x_cvbdyshape
    t_cvbdyshape_full=x_cvbdyshape

    call incref(x_cvshape_full)
    call incref(t_cvshape_full)
    call incref(x_cvbdyshape_full)
    call incref(t_cvbdyshape_full)
  end if

  call allocate(cfl_no, tfield%mesh, name="DummyCFLNumber", &
                field_type=FIELD_TYPE_CONSTANT)
  call zero(cfl_no)
  fdcflno = .false.

  ! get the mesh sparsity for the matrices
  if(tfield_options%diffusionscheme==CV_DIFFUSION_BASSIREBAY) then
    ! extend the sparsity
    mesh_sparsity=make_sparsity_transpose(tfield%mesh, diffusivity%mesh, "AdvectionDiffusionSparsity")

    mesh_sparsity_1=make_sparsity(tfield%mesh, tfield%mesh, "AdvectionSparsity")

    if(.not.(tfield%mesh==diffusivity%mesh)) then
      if(tfield%mesh%shape%degree>1) then
        FLAbort("To have a different diffusivity mesh the field must be at most P1")
      elseif(diffusivity%mesh%shape%degree>tfield%mesh%shape%degree) then
        FLAbort("The diffusivity mesh must be of a lower degree than the field")
      end if

      grad_m_t_sparsity=make_sparsity(tfield%mesh, diffusivity%mesh, "GradientTransposedSparsity")
    else
      grad_m_t_sparsity=mesh_sparsity_1
      call incref(grad_m_t_sparsity)
    end if

  else

    mesh_sparsity=make_sparsity(tfield%mesh, tfield%mesh, "AdvectionDiffusionSparsity")

    grad_m_t_sparsity=mesh_sparsity
    call incref(grad_m_t_sparsity)

    mesh_sparsity_1=mesh_sparsity
    call incref(mesh_sparsity_1)
  end if

  ! allocate the advection matrix
  call allocate(A_m, mesh_sparsity, name="AdvectionMatrix")
  call zero(A_m)

  call allocate(D_m, sparsity=mesh_sparsity, name="AuxiliaryMatrix")
  call zero(D_m)

  call allocate(inverse_mass, sparsity=mesh_sparsity_1, name="InverseMassMatrix")
  call zero(inverse_mass)

  call allocate(mass, sparsity=mesh_sparsity_1, name="MassMatrix")
  call zero(mass)

  call allocate(diff_rhs, tfield%mesh, name="DiffusionRHS")
  call zero(diff_rhs)

  call allocate(rhs, tfield%mesh, "RHS")
  call zero(rhs)

  if(tfield_options%diffusionscheme==CV_DIFFUSION_BASSIREBAY) then
    call allocate(q_lumpedmass, diffusivity%mesh, name="DiffusivityLumpedMass")
    call compute_lumped_mass(x, q_lumpedmass)
  else
    call allocate(q_lumpedmass, diffusivity%mesh, name="DiffusivityLumpedMass", field_type=FIELD_TYPE_CONSTANT)
    call zero(q_lumpedmass)
  end if

  call allocate(lumpedmass, tfield%mesh, name="FieldLumpedMass")
  call compute_lumped_mass(x, lumpedmass)
  call addto_diag(mass, lumpedmass)

  lumpedmass%val = 1./lumpedmass%val
  call addto_diag(inverse_mass, lumpedmass)

  getmat = .true.

  call assemble_advectiondiffusion_m_cv(A_m, rhs, &
                              tfield, old_tfield, tfield_options, &
                              tdensity, old_tdensity, tdensity_options, &
                              cvfaces, x_cvshape, x_cvbdyshape, &
                              u_cvshape, u_cvbdyshape, t_cvshape, &
                              states, advu, x, cfl_no, fdcflno, &
                              getmat, dt, &
                              mesh_sparsity_1, &
                              diffusion=diffusion, diffusivity=diffusivity, q_lumpedmass=q_lumpedmass, &
                              D_m=D_m, diff_rhs=diff_rhs, grad_m_t_sparsity=grad_m_t_sparsity, &
                              x_cvshape_full=x_cvshape_full, x_cvbdyshape_full=x_cvbdyshape_full, &
                              t_cvshape_full=t_cvshape_full, t_cvbdyshape_full=t_cvbdyshape_full)

  call apply_dirichlet_conditions(mass, rhs, tfield, dt)

  call mmwrite(trim(simname)//"_diffusion.mm", D_m)
  call mmwrite(trim(simname)//"_mass.mm", mass)
  call mmwrite(trim(simname)//"_inversemass.mm", inverse_mass)

  call deallocate(states(1))
  call deallocate(rhs)
  call deallocate(mesh_sparsity)
  call deallocate(mesh_sparsity_1)
  call deallocate(grad_m_t_sparsity)
  call deallocate(A_m)
  call deallocate(D_m)
  call deallocate(cfl_no)
  call deallocate(diff_rhs)
  call deallocate(q_lumpedmass)
  call deallocate(x_cvshape)
  call deallocate(u_cvshape)
  call deallocate(t_cvshape)
  call deallocate(x_cvbdyshape)
  call deallocate(u_cvbdyshape)
  call deallocate(x_cvshape_full)
  call deallocate(t_cvshape_full)
  call deallocate(x_cvbdyshape_full)
  call deallocate(t_cvbdyshape_full)
  call deallocate(advu)
  call deallocate(dummydensity)
  call deallocate(cvfaces)
  call deallocate(inverse_mass)
  call deallocate(mass)
  call deallocate(lumpedmass)

  call print_references(0)

#ifdef HAVE_MPI
  call MPI_Finalize(ierr)
  assert(ierr == MPI_SUCCESS)
#endif

contains

  subroutine read_command_line(filename)
    character(len=*), intent(out) :: filename
    integer :: status

    call get_command_argument(1, value=filename, status=status)

    select case(status)
    case(1:)
       call usage
       stop
    case(:-1)
       write(0,*) "Warning: truncating filename"
    end select

  end subroutine read_command_line

  subroutine usage
    write(0,*) "usage:"
    write(0,*) "cv_diffusion_matrix flml-file"
    write(0,*) "Dumps the cv diffusion matrix."
  end subroutine usage

end program cv_diffusion_matrix
