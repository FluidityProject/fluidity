      Module mba3d_auxSP
C
      use mba3d_error
      use mba3d_findSE
C
      contains
C
C ================================================================
      Subroutine makSP(iP, IEP, IPE, IEE, MaxS, lE, iEs)
C ================================================================
C Remark: the 3rd column of IPE is overloaded in makSE. If it is
C         needed, the absolute value should be used.
c 
c         We realize the same idea here by overloading the 3rd 
c         column of IPE. 
C ================================================================
c group (M)
      Integer IEP(*)
      Integer IPE(5, *), IEE(4, *)

C group (S)
      Integer iEs(*)

C group (Local variables)
      Logical repeat

C ================================================================
      lE = 1
      iEs(lE) = IEP(iP)

      iE = iEs(1)
      IPE(3, iE) = -IPE(3, iE)

      n2 = 0
 1    repeat = .FALSE.

      n1 = n2 + 1
      n2 = lE
      Do 4 n = n1, n2
         If(lE.GE.MaxS - 4) Goto 1000

         Do 2 i = 1, 4
            iE = IEE(i, iEs(n))
            If(iE.EQ.0) Goto 2
            If(IPE(3, iE).LT.0) Goto 2

            If(iP.EQ.IPE(1, iE) .OR.
     &         iP.EQ.IPE(2, iE) .OR.
     &         iP.EQ.IPE(3, iE) .OR.
     &         iP.EQ.IPE(4, iE)) Then
c              Do k = lE, 1, -1
c                 If(iE.EQ.iEs(k)) Goto 2
c              End do

               repeat = .TRUE.
               lE = lE + 1

               iEs(lE) = iE
               IPE(3, iE) = -IPE(3, iE)
            End if
 2       Continue
 4    Continue

      If(repeat) Goto 1

c ... restoring the overloaded values
      Do k = 1, lE
         iE = iEs(k)
         IPE(3, iE) = -IPE(3, iE)
      End do

      Return
 1000 Call errMes(1007, 'makSP', 'local variable MaxS is small')
      End Subroutine makSP



C ================================================================
      Subroutine chkSPf(NM, iP1, iOPERAT, ICP, IEP, IPE, IEE, lPf, iPf)
C ================================================================
C  Routine collects boundary points within NM mesh steps of iP1.
C ================================================================
      include 'makS.fd'
      include 'color.fd'
      include 'operat.fd'
C ================================================================
c group (M)
      Integer ICP(*), IEP(*)
      Integer IPE(5, *), IEE(4, *)

C group (S)
      Integer iPf(*)

C group (Local variables)
      Integer iEs(MaxS)

C ================================================================
      If(NM.GT.2) Call errMes(6001, 'chkSPf', 'system error')

      lPf = 0
      Do nStep = 1, NM
         If(nStep.EQ.1) Then
            m2 = 1
            iPf(1) = iP1
         Else
            m2 = lPf
         End if

         Do m = 1, m2
            Call makSP(iPf(m), IEP, IPE, IEE, MaxS, lE, iEs)
            Do n = 1, lE
               iE = iEs(n)
               Do 10 i = 1, 4
                  iPt = IPE(i, iE)
                  If(iOPERAT.EQ.iDELET .AND. iPt.EQ.iP1) Goto 10

                  If(ifXnode(ICP(iPt), jBnode)) Then
                     Call findSE(lPf, iPf, iPt, nPt)
                     If(nPt.NE.0) Goto 10

                     lPf = lPf + 1
                     If(lPf.GT.MaxS) Goto 1000
                     iPf(lPf) = iPt
                  End if
 10            Continue
            End do
         End do
      End do

      Return
 1000 Call errMes(1007, 'chkSPf', 'local variable MaxS is small')
      End Subroutine chkSPf



C ================================================================
      Subroutine info
C ================================================================
C Routines prints ANI's logo. 
C ================================================================
      Integer iDat(14)
      DATA iDat/83,  116, 111, 110, 101, 32, 70,
     &          108, 111, 119, 101, 114, 32, 33/

      Write(*,'(14A1)') (char(iDat(i)), i = 1, 14)
      Return
      End Subroutine info



C ================================================================
      Subroutine chkSPb(
C ================================================================
     &           NM, lD, iDs, lN, iNs, iOPERAT,
     &           ICP, IEP, IPE, IEE, lPf, iPf, iPbad, flag)
C ================================================================
C  Routine checks that each boundary point in array lPf can be
C  connected with an interior point (having color jInode) using
C  at most NM mesh edges. flag = T, when there is a point which
C  is surrounded by boundary points.
C
C  The first bad point is returned in iPbad.
C ================================================================
      include 'makS.fd'
      include 'color.fd'
      include 'operat.fd'
C ================================================================
c group (M)
      Integer ICP(*), IEP(*)
      Integer IPE(5, *), IEE(4, *)

C group (S)
      Integer iDs(2, *), iNs(2, *), iPf(*)
      Logical flag

C group (Local variables)
      Integer iPb(MaxS), iEs(MaxS)

C ================================================================
      If(NM.GT.2) Call errMes(6001, 'chkSPf', 'system error')

      flag = .TRUE.
      Do 10 l = 1, lPf
         lPb = 1
         iPb(1) = iPf(l)

         Do nStep = 1, NM
            m2 = lPb

            Do m = 1, m2
               iPc = iPb(m)
               Call makSP(iPc, IEP, IPE, IEE, MaxS, lE, iEs)
               Do n = 1, lE
                  iE = iEs(n)
                  Do 5 i = 1, 4
                     iPt = IPE(i, iE)

                     Do k = 1, lD
                        If(iPc.EQ.iDs(1, k) .AND.
     &                     iPt.EQ.iDs(2, k)) Goto 5
                        If(iPc.EQ.iDs(2, k) .AND.
     &                     iPt.EQ.iDs(1, k)) Goto 5
                     End do

                     If(ifXnode(ICP(iPt), jInode)) Goto 10 

                     Call findSE(lPb, iPb, iPt, nPt)
                     If(nPt.NE.0) Goto 5

                     lPb = lPb + 1
                     If(lPb.GT.MaxS) Goto 1000
                     iPb(lPb) = iPt

                     Do k = 1, lN
                        Do j = 1, 2
                           If(iPc.EQ.iNs(j, k)) Then
                              iPt = iNs(3 - j, k)
                              If(ifXnode(ICP(iPt), jInode)) Goto 10

                              If(nStep.LT.NM) Then
                                 Call findSE(lPb, iPb, iPt, nPt)
                                 If(nPt.NE.0) Goto 5

                                 lPb = lPb + 1
                                 If(lPb.GT.MaxS) Goto 1000
                                 iPb(lPb) = iPt
                              End if

c  ...  the next line restricted the search for inner nodes and was removed
c                              Goto 5
                           End if
                        End do
                     End do
 5                Continue
               End do
            End do
         End do

         iPbad = iPb(1)
         Return
 10   Continue

      flag = .FALSE.
      Return
 1000 Call errMes(1007, 'chkSP', 'local variable MaxS is small')
      End Subroutine chkSPb


C ================================================================
      Real*8 Function calNorm(xyz)
C ================================================================
C Routines computes 2-norm of vector xyz.
C ================================================================
      Real*8 xyz(3)

      calNorm = dsqrt(xyz(1) ** 2 + xyz(2) ** 2 + xyz(3) ** 2)

      Return
      End Function calNorm


C ================================================================
      Real*8 Function sqrNorm(xyz)
C ================================================================
C Routines computes square of 2-norm of vector xyz.
C ================================================================
      Real*8 xyz(3)

      sqrNorm = xyz(1) ** 2 + xyz(2) ** 2 + xyz(3) ** 2

      Return
      End Function sqrNorm


C ================================================================
C  All operations with colors are binary operations. 
C  The same operations go with status.
C ================================================================


C ================================================================
C     Logical Function ifXstatus(status, iXstatus)
      Logical Function ifXnode(clr, iXnode)
C ================================================================
      Integer clr, iXnode

      ifXnode = IAND(clr, iXnode) .EQ. iXnode

      Return
      End Function ifXnode



C ================================================================
C     Subroutine addXstatus(status, iXstatus)
      Subroutine addXnode(clr, iXnode)
C ================================================================
C Binary operation +.
C ================================================================
      Integer clr, iXnode

      clr = IOR(clr, iXnode)

      Return
      End Subroutine addXnode



C ================================================================
C     Subroutine delXstatus(status, iXstatus)
      Subroutine delXnode(clr, iXnode)
C ================================================================
      Integer clr, iXnode

      clr = clr - IAND(clr, iXnode)

      Return
      End Subroutine delXnode



C ================================================================
      Integer Function minClr(clr1, clr2)
C ================================================================
C  The function returns common color for both clr1 and clr2.
C ================================================================
      Integer clr1, clr2

      minClr = IAND(clr1, clr2)

      Return
      End Function minClr



C ================================================================
      Integer Function maxClr(clr1, clr2)
C ================================================================
C  The function returns minimal color containing both clr1 and clr2.
C ================================================================
      Integer clr1, clr2

      maxClr = IOR(clr1, clr2)

      Return
      End Function maxClr



C ================================================================
      Subroutine setStatus(flagAuto, status, iPrint)
C ================================================================
      include 'status.fd'
C ================================================================
C Routine adds additional properties to the variable status
C ================================================================
      Logical  flagAuto
      Integer  status, iPrint

C ================================================================
      status = max(0, status)

c ... remove obsolete and not-implemented input features 
      If(status.GT.0) Then
         If(ifXnode(status, ANITangledMesh)) Then
            Call delXnode(status, ANITangledMesh)
            If(iPrint.GE.2) Write(*, 5001) 
         End if
      End if

      
c ... inform the user about requested features  
      If(iPrint.GE.2) Then
         If(ifXnode(status, ANIForbidBoundaryElements)) Write(*, 5003) 
         If(ifXnode(status, ANIFixBoundaryFaces))       Write(*, 5004)
         If(ifXnode(status, ANIFixSurfacePoints))       Write(*, 5007)

         If(ifXnode(status, ANIMultiConnectedGeometry)) 
     &      Write(*, 5002) '[user]'
         If(ifXnode(status, ANIUse2ArmRule))      
     &      Write(*, 5005) '[user]'
         If(ifXnode(status, ANIDeleteTemporaryFaces))
     &      Write(*, 5006) '[user]'
      End if


c ... set up default features
      If(flagAuto) Then
         If(.NOT.ifXnode(status, ANIUse2ArmRule)) Then
            If(iPrint.GE.2) Write(*, 5005) '[system]'

            Call addXnode(status, ANIUse2ArmRule)
         End if

         If(.NOT.ifXnode(status, ANIDeleteTemporaryFaces)) Then
            If(iPrint.GE.2) Write(*, 5006) '[system]'

            Call addXnode(status, ANIDeleteTemporaryFaces)
         End if
      End if

      If(iPrint.GE.2) Write(*,*) 

      Return

 5001 Format('status.fd: -512  [ANITangledMesh]            [research]')
 5002 Format('status.fd: +64   [ANIMultiConnectedGeometry] ', A)
 5003 Format('status.fd: +1    [ANIForbidBoundaryElements] [user]') 
 5004 Format('status.fd: +4    [ANIFixBoundaryFaces]       [user]')
 5005 Format('status.fd: +2    [ANIUse2ArmRule]            ', A)
 5006 Format('status.fd: +8    [ANIDeleteTemporaryFaces]   ', A)
 5007 Format('status.fd: +16   [ANIFixSurfacePoints]       [user]')

      Return
      End Subroutine setStatus
C
      End Module mba3d_auxSP
