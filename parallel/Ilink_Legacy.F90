! Not in a module as called from C++

! This function remains to be called from c (gets mapped to ilink2_fc->ilink2)
! until flcomms has been sufficiently generalised.
function ilink2_legacy(iloc, jloc) result(ilink2)
  ! Returns the local node number in a quadratic tet, of the node on the 
  ! edge between local nodes iloc and jloc, or the vertex if iloc==jloc.
  ! iloc and jloc are in the local linear tet numbering. The result
  ! is in the 'one true element numbering' for quadratic tets.
  integer ilink2
  integer, intent(in):: iloc, jloc

  select case (iloc)
  case (1)
     if (jloc==1) ilink2=1
     if (jloc==2) ilink2=2
     if (jloc==3) ilink2=4
     if (jloc==4) ilink2=7
  case (2)
     if (jloc==1) ilink2=2
     if (jloc==2) ilink2=3
     if (jloc==3) ilink2=5
     if (jloc==4) ilink2=8
  case (3)
     if (jloc==1) ilink2=4
     if (jloc==2) ilink2=5
     if (jloc==3) ilink2=6
     if (jloc==4) ilink2=9
  case (4)
     if (jloc==1) ilink2=7
     if (jloc==2) ilink2=8
     if (jloc==3) ilink2=9
     if (jloc==4) ilink2=10
  end select

end function ilink2_legacy
