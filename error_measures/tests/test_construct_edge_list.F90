subroutine test_construct_edge_list

  use elements
  use state_module
  use linked_lists
  use sparse_tools
  use fields
  use gradation_metric
  use unittest_tools
  use adjacency_lists
  use vtk_interfaces
  implicit none

  type(state_type) :: state
  type(mesh_type), pointer :: mesh
  type(csr_sparsity), pointer :: nn_sparsity
  type(csr_matrix) :: nnlist
  type(elist) :: edgelist
  logical :: fail
  integer :: correct_size

  call vtk_read_state("data/cube.vtu", state)
  mesh => extract_mesh(state, "Mesh")

  nn_sparsity => extract_nnlist(mesh)
  call allocate(nnlist, nn_sparsity, type=CSR_INTEGER)
  nnlist%ival = -1

  call construct_edge_list(mesh, nnlist, edgelist)

  correct_size = (size(nnlist%ival)) / 2              ! the number of nonzero
                                                      ! entries in the upper
                                                      ! triangle of the matrix
  fail = .false.
  if (edgelist%length /= correct_size) then
    write(0,*) "edgelist%length == ", edgelist%length
    write(0,*) "correct_size == ", correct_size
    fail = .true.
  end if

  call report_test("[edgelist size]", fail, .false., "Edgelist should be the &
  & correct size.")

  fail = .false.
  if (any(nnlist%ival == -1)) fail = .true.
  call report_test("[edgelist complete]", fail, .false., "Edgelist should &
  & have every edge in it.")

  deallocate(nnlist%ival)
  call flush_list(edgelist)
  call deallocate(state)

end subroutine test_construct_edge_list
