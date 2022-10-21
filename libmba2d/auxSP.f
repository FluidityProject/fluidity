C ================================================================
C  All operations with colors are binary operations. 
C  The same operations go with the input parameter status.
C ================================================================


C ================================================================
C     Logical Function ifXstatus(status, iXstatus)
      Logical Function ifXnode(clr, iXnode)
C ================================================================
      Integer clr, iXnode

      ifXnode = IAND(clr, iXnode) .EQ. iXnode

      Return
      End



C ================================================================
C     Subroutine addXstatus(status, iXstatus)
      Subroutine addXnode(clr, iXnode)
C ================================================================
      Integer clr, iXnode

      clr = IOR(clr, iXnode)

      Return
      End



C ================================================================
C     Subroutine delXstatus(status, iXstatus)
      Subroutine delXnode(clr, iXnode)
C ================================================================
      Integer clr, iXnode

      clr = clr - IAND(clr, iXnode)

      Return
      End



C ================================================================
      Integer Function minClr(clr1, clr2)
C ================================================================
C  The function returns common color for both clr1 and clr2.
C ================================================================
      Integer clr1, clr2

      minClr = IAND(clr1, clr2)

      Return
      End



C ================================================================
      Integer Function maxClr(clr1, clr2)
C ================================================================
C Routine returns minimal color containing both clr1 and clr2.
C ================================================================
      Integer clr1, clr2

      maxClr = IOR(clr1, clr2)

      Return
      End



C ================================================================
      Subroutine setStatus(flagAuto, status, iPrint)
C ================================================================
      include 'status.fd'
C ================================================================
C Routine adds additional properties to the variable status
C ================================================================
      Logical  flagAuto
      Integer  status, iPrint

      Logical  ifXnode

C ================================================================
      status = max(0, status)

c ... remove obsolete and not-implemented input features 
      If(status.GT.0) Then
         If(ifXnode(status, ANISmoothMesh)) Then
            Call delXnode(status, ANISmoothMesh)
            If(iPrint.GE.1) Write(*, 5001) 
         End if

         If(ifXnode(status, ANIMultiConnectedGeometry)) Then
           Call delXnode(status, ANIMultiConnectedGeometry)
            If(iPrint.GE.1) Write(*, 5002) 
         End if
      End if

      
c ... inform the user about requested features  
      If(iPrint.GE.1) Then
         If(ifXnode(status, ANIForbidBoundaryElements)) Write(*, 5003) 
         If(ifXnode(status, ANIFixBoundaryEdges))       Write(*, 5004)
         If(ifXnode(status, ANIFixBoundaryPoints))      Write(*, 5007)

         If(ifXnode(status, ANIUse2ArmRule))      
     &      Write(*, 5005) '[user]'
         If(ifXnode(status, ANIDeleteTemporaryEdges))
     &      Write(*, 5006) '[user]'
         If(ifXnode(status, ANIUntangleMesh))
     &      Write(*, 5008) '[user]'
      End if


c ... set up default features
      If(flagAuto) Then
         If(.NOT.ifXnode(status, ANIUse2ArmRule)) Then
            If(iPrint.GE.2) Write(*, 5005) '[system]'

            Call addXnode(status, ANIUse2ArmRule)
         End if

         If(.NOT.ifXnode(status, ANIDeleteTemporaryEdges)) Then
            If(iPrint.GE.2) Write(*, 5006) '[system]'

            Call addXnode(status, ANIDeleteTemporaryEdges)
         End if
      End if

      If(iPrint.GE.1) Write(*,*) 

      Return

 5001 Format('status.fd: -256  [ANISmoothMesh]             [obsolete]')
 5002 Format('status.fd: -64   [ANIMultiConnectedGeometry] [not used]')
 5003 Format('status.fd: +1    [ANIForbidBoundaryElements] [user]') 
 5004 Format('status.fd: +4    [ANIFixBoundaryEdges]       [user]')
 5005 Format('status.fd: +2    [ANIUse2ArmRule]            ', A)
 5006 Format('status.fd: +8    [ANIDeleteTemporaryEdges]   ', A)
 5007 Format('status.fd: +16   [ANIFixBoundaryPoints]      [user]')
 5008 Format('status.fd: +512  [ANIUntangleMesh]           ', A)

      Return
      End


