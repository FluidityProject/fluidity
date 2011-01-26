#include "fdebug.h"

module zoltan_callbacks

  use zoltan
  use zoltan_global_variables

  ! Needed for zoltan_cb_owned_node_count
  use halos_ownership
  use metric_tools  

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
 
end module zoltan_callbacks
