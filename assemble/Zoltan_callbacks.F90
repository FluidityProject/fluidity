#include "fdebug.h"

module zoltan_callbacks

  use zoltan
  use zoltan_global_variables

  ! Needed for zoltan_cb_owned_node_count
  use halos, only: halo_nowned_nodes, get_owned_nodes, halo_universal_number
  use metric_tools

  ! Needed for zoltan_cb_get_owned_nodes
  ! - added get_owned_nodes, halo_universal_number to use halos
  use sparse_tools, only: row_length

  implicit none
  
  public
  
contains

  function zoltan_cb_owned_node_count(data, ierr) result(count)
    integer(zoltan_int) :: count
    integer(zoltan_int), dimension(*) :: data ! not used
    integer(zoltan_int), intent(out) :: ierr
    
    ewrite(1,*) "In zoltan_cb_owned_node_count"

    count = halo_nowned_nodes(zoltan_global_zz_halo)
    if (have_option("/mesh_adaptivity/hr_adaptivity/zoltan_options/zoltan_debug")) then
       ewrite(1,*) "zoltan_cb_owned_node_count found: ", count, " nodes"
    end if
    ierr = ZOLTAN_OK
  end function zoltan_cb_owned_node_count


  subroutine zoltan_cb_get_owned_nodes(data, num_gid_entries, num_lid_entries, global_ids, local_ids, wgt_dim, obj_wgts, ierr)
    integer(zoltan_int), dimension(*), intent(in) :: data ! not used
    integer(zoltan_int), intent(in) :: num_gid_entries, num_lid_entries 
    integer(zoltan_int), intent(out), dimension(*) :: global_ids 
    integer(zoltan_int), intent(out), dimension(*) :: local_ids 
    integer(zoltan_int), intent(in) :: wgt_dim 
    real(zoltan_float), intent(out), dimension(*) :: obj_wgts 
    integer(zoltan_int), intent(out) :: ierr
    
    integer :: count, i
    real(zoltan_float) :: max_obj_wgt
    
    ewrite(1,*) "In zoltan_cb_get_owned_nodes"
    
    assert(num_gid_entries == 1)
    assert(num_lid_entries == 1)
    assert(wgt_dim == 1)
    
    count = halo_nowned_nodes(zoltan_global_zz_halo)
    
    call get_owned_nodes(zoltan_global_zz_halo, local_ids(1:count))
    global_ids(1:count) = halo_universal_number(zoltan_global_zz_halo, local_ids(1:count))
    
    if (have_option("/mesh_adaptivity/hr_adaptivity/zoltan_options/zoltan_debug")) then
       ewrite(1,*) "zoltan_cb_get_owned nodes found local_ids: ", local_ids(1:count)
       ewrite(1,*) "zoltan_cb_get_owned nodes found global_ids: ", global_ids(1:count)
    end if
    
    if(zoltan_global_migrate_extruded_mesh) then
       ! weight the nodes according to the number of nodes in the column beneath it
       max_obj_wgt = 1.0
       do i=1,count
          obj_wgts(i) = float(row_length(zoltan_global_columns_sparsity, i))
          max_obj_wgt = max(max_obj_wgt, obj_wgts(i))
       end do
       ! normalise according to the most nodes in a column
       do i=1,count
          obj_wgts(i) = obj_wgts(i)/max_obj_wgt
       end do
    else
       do i=1,count
          obj_wgts(i) = 1.0
       end do
    end if
    
    ierr = ZOLTAN_OK
  end subroutine zoltan_cb_get_owned_nodes
 
end module zoltan_callbacks
