      Module mba3d_auxSR
C
      use mba3d_auxSE
      use mba3d_auxSP
C
      contains
C
C ================================================================
      Subroutine makSR(iPa, iPb, lE, iEs, IPEs, lR, iRs, flag)
C ================================================================
C The ordered sequence of edges "orthogonal" to the edge {iPa, iPb}
C is computed.
C
C Remark: The superlement is not modified.
C ================================================================
c group (S)
      Integer iEs(*)
      Integer IPEs(5, *), iRs(3, *)

      Logical flag

C group (Local variables)

C ================================================================
      flag = .TRUE.
      lR = 0
      Do n = 1, lE
         iP1 = IPEs(1, n)
         iP2 = IPEs(2, n)
         iP3 = IPEs(3, n)
         iP4 = IPEs(4, n)

         If(check14(iPa, iP1, iP2, iP3, iP4) .AND.
     &      check14(iPb, iP1, iP2, iP3, iP4)) Then
            lR = lR + 1
            iRs(1, lR) = n

            icnt = 1
            Do i = 1, 4
               iPt = IPEs(i, n)
               If(iPt.NE.iPa .AND. iPt.NE.iPb) Then
                  icnt = icnt + 1
                  iRs(icnt, lR) = iPt
               End if
            End do
         End if
      End do


c ... reodering of the edges (2 iterations are allowed)
      itr = 0
 10   Do n = 1, lR - 1
         iP1 = iRs(2, n)
         iP2 = iRs(3, n)

         Do m = n + 1, lR
            iP3 = iRs(2, m)
            iP4 = iRs(3, m)

            m1 = m
            If(iP2.EQ.iP3) Then
               Goto 20
            Else if(iP2.EQ.iP4) Then
               iRs(2, m) = iP4
               iRs(3, m) = iP3
               Goto 20
            End if
         End do

         Do i = 1, 3
            iwk = iRs(i, 1)
            iRs(i, 1) = iRs(i, n)
            iRs(i, n) = iwk
         End do
         iRs(2, 1) = iP2
         iRs(3, 1) = iP1

         itr = itr + 1
         If(itr.EQ.3) Then
c  ...  there is no connected path
            flag = .FALSE.
            Goto 1000
         End if
         Goto 10


 20      Do i = 1, 3
            iwk = iRs(i, n + 1)
            iRs(i, n + 1) = iRs(i, m1)
            iRs(i, m1) = iwk
         End do
      End do
 1000 Return
      End Subroutine makSR


C ================================================================
      Subroutine clrSR(iPa, iPb, ICP, IPF, IFE,
     &                 lF, iFs, lE, iEs, ICRab)
C ================================================================
C The color of edge {iPa, iPb} is computed. 
C
C Remark: The superlement is not modified.
C ================================================================
      include 'makS.fd'
      include 'color.fd'
C ================================================================
      Integer  ICP(*), IPF(4, *), IFE(4, *), iFs(*), iEs(*)

C group (Local variables)
      Integer ip(4), icnt(MaxS)
      Logical flag

C ================================================================
      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 1

      flag = .FALSE.

      Do n = 1, lF
         icnt(n) = -1
      End do

      iCRab = 0
      Do 10 n = 1, lF
         iFt = iFs(n)
         If(iFt.LE.0) Goto 10

         Do i1 = 1, 3
            i2 = ip(i1 + 1)

            iP1 = IPF(i1, iFt)
            iP2 = IPF(i2, iFt)
            If(iPa.EQ.iP1 .AND. iPb.EQ.iP2 .OR.
     &         iPb.EQ.iP1 .AND. iPa.EQ.iP2) Then
               iCLRf = IPF(4, iFt)
               icnt(n) = 0

               If(iCRab.EQ.0) Then
                  iFp = iFt

                  iCRab = iCLRf
               Else If(iCRab.NE.iCLRf) Then
                  iCRab = minClr(ICP(iPa), ICP(iPb))
c  ...  a temporary fix for V-V edges
                  Call delXnode(iCRab, jVnode)

                  Call addXnode(iCRab, jRnode)
                  Goto 1000
               End if
            End if
         End do
 10   Continue


      If(iCRab.NE.0) Then
         iCRab = jSnode

c  ...  checking for edges on boundary surfaces
         ic1 = 0
         Do k = 1, lE
            iEt = iEs(k)
            If(iEt.GT.0) Then
               Do i = 1, 4
                  If(IFE(i, iEt).EQ.iFp) ic1 = ic1 + 1
                  If(ic1.EQ.2) Goto 1000
               End do
            End if
         End do

         Call addXnode(iCRab, jBnode)
      Else
         ICRab = jInode
      End if
 1000 Return
      End Subroutine clrSR



C ================================================================
      Real*8 Function angle2Edges( XYA, XYB, XYC )
C ================================================================
C Cosine of angle between two edges {A,B} and {A,C}.
C ================================================================
      Real*8  XYA(3), XYB(3), XYC(3)
      Real*8  XYN(3), XYM(3)
      Real*8  d1, d2

      Do i = 1, 3
         XYN(i) = XYB(i) - XYA(i)
         XYM(i) = XYC(i) - XYA(i)
      End do

      d1 = XYN(1) ** 2 + XYN(2) ** 2 + XYN(3) ** 2
      d2 = XYM(1) ** 2 + XYM(2) ** 2 + XYM(3) ** 2

      angle2Edges = (XYN(1) * XYM(1) +
     &               XYN(2) * XYM(2) +
     &               XYN(3) * XYM(3)) / dsqrt(d1 * d2)

      Return
      End Function angle2edges



C ================================================================
      Real*8 Function calEdge(xy1, xy2)
C ================================================================
      Real*8 xy1(3), xy2(3)

      calEdge = dsqrt((xy1(1) - xy2(1)) ** 2 +
     &                (xy1(2) - xy2(2)) ** 2 +
     &                (xy1(3) - xy2(3)) ** 2)
      Return
      End Function calEdge



C ================================================================
      Real*8 Function sqrEdge(xy1, xy2)
C ================================================================
      Real*8 xy1(3), xy2(3)

      sqrEdge = (xy1(1) - xy2(1)) ** 2 +
     &          (xy1(2) - xy2(2)) ** 2 +
     &          (xy1(3) - xy2(3)) ** 2
      Return
      End Function sqrEdge
C
      End Module mba3d_auxSR
