subroutine test_bubble_tools

  use fields
  use state_module
  use sparse_matrices_fields
  use vtk_interfaces
  use unittest_tools
  use quadrature
  use bubble_tools
  use spud
  use populate_state_module

  type(vector_field), pointer :: X
  logical :: fail = .false., warn = .false.
  type(quadrature_type) :: quad
  type(element_type) :: shape
  integer :: i
  real, dimension(7) :: N_Vals
  type(state_type), dimension(:), pointer :: states
  type(mesh_type), pointer :: cg_mesh, dg_mesh
  type(scalar_field) :: p2bubble_field, p1dg_field2
  type(scalar_field), pointer :: p1dg_field
  type(csr_matrix), pointer :: cgdg_projection

  quad = make_quadrature(vertices=3,dim=2,degree=8,family=FAMILY_COOLS)
  shape = make_element_shape(vertices=3,dim=2,degree=2,quad=quad,&
       &type=ELEMENT_BUBBLE)
  fail = shape%loc.ne.7
  call report_test("[P2b DOF count]", fail, .false., "P2b should have 7 DOFS") 

  N_vals = eval_shape(shape, (/1.0/3.0,1.0/3.0,1.0/3.0/))
  fail = maxval(abs(N_vals-(/ -0.11111111111111110, 0.44444444444444442,&
       &-0.11111111111111110, 0.44444444444444442, 0.44444444444444442,&
       &-0.11111111111111110, 3.70370370370370350E-002/)))>1.0e-8
  call report_test("[P2b barycentre]", fail,&
       &.false., "Wrong values at barycentre")

  call nodalise_bubble_basis(shape)
  fail = any(abs(sum(shape%n,1)-1.0)>1.0e-8)
  call report_test("[Partition of unity]", fail,&
       &.false., "Should sum to 1.")  

  call load_options("data/blob.flml")
  call populate_state(states)

  ! set up projection
  cg_mesh => extract_mesh(states(1),"P2BubbleMesh")
  dg_mesh => extract_mesh(states(1),"DGMesh")
  call setup_Cg_Dg_projection(states(1),Cg_mesh,Dg_mesh)

  ! project p1dg field into a p2bubble field
  call allocate(P2bubble_field,cg_mesh,"P2bubble_field")
  p1dg_field => extract_scalar_field(states(1),name="P1DGField")
  cgdg_projection => extract_csr_matrix(states(1),&
       &trim(Cg_mesh%name)//trim(Dg_mesh%name)//"Projection")
  call mult(P2bubble_field,cgdg_projection,p1dg_field)

  ! try to undo the projection
  call allocate(p1dg_field2,dg_mesh,"P1DGField2")
  call cg_dg_projection(states(1),P2bubble_field,p1dg_field2,1.0e-6)
end subroutine test_bubble_tools
