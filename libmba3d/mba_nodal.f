      Module mba3d_mba_nodal
C
      use mba3d_ani2
      use mba3d_check
      use mba3d_forlibmba
C
      contains
C
C ================================================================
      Subroutine mbaNodal(
C ================================================================
c group (M)
     &      nP, MaxP, nF, MaxF, nE, MaxE,
     &      XYP, IPF, IPE, lbF, lbE,
     &      nEStar, 
c group (Dev)
     &      nPv, nFv, nEv, IPV, IFV, IEV, 
     &      flagAuto, status,
c group (Q)
     &      MaxSkipE, MaxQItr,
     &      Metric, Quality, rQuality, 
c group (W)
     &      MaxWr, MaxWi, rW, iW,
     &      iPrint, iERR)
C ================================================================
      include 'status.fd'
      include 'magic.fd'
      include 'output.fd'
C ================================================================
C  VARIABLES & PARAMETER:
C     All input arguments are divided into two groups. The first
C     group (BASIC) containes the main parameters. The second
C     group (AUX) contains auxiliary parameters which can be
C     omitted.
C
C  ***BASIC 1: nP to IPF
C     nP, MaxP  - actual and maximal number of points (P)
C     nF, MaxF  - actual and maximal number of boundary and inner faces (F)
C     nE, MaxE  - actual and maximal number of elements (E)
C
C     XYP(3, MaxP) - Cartesian coordinates of mesh points
C     IPE(4, MaxE) - connectivity list of element
C     IPF(3, MaxF) - connectivity list of faces
C
C  ***AUX 1: lbF to lbE. 
C     lbF(MaxF) - boundary identificator
C                   (a positive number strictly less than 1400=MaxS-100,
C                    colors from MaxS-99 to MaxS are reserved)
C                   (example: unit cube has 6 boundaries which
C                    may have or not different identificators)
C
C     lbE(MaxE) - element indentificator (a positive number)
C
C  ***BASIC 2: nEStar
C     nEstar    - target number of elements in the final mesh
C
C  ***AUX 2: nPv to MaxSkipE
C     nPv       - number of fixed points     (can be zero)
C     nFv       - number of fixed triangles  (can be zero)
C     nEv       - number of fixed tetrahedra (can be zero)
C
C     IPV(nPv)  - list of fixed points
C     IFV(nFv)  - list of fixed triangles
C     IEV(nEv)  - list of fixed tetrahedra
C
C     flagAuto  - flag controling the mesh generation:
C                 TRUE  - automatic recovering of missing mesh elements
C                 FALSE - rigorous checking of input data
C
C     status    - advanced control of mesh generation:
C                 0 or negative - no additional control
C                 positive      - enforce some of the mesh properties
C                    The detailed description of available properties
C                    is in file status.fd. Variable status is equal
C                    to the sum of positive numbers corresponding
C                    to the desired properties. Here is a short list
C                    of user-controled properties:
C
C                     1 - final mesh will not have boundary elements
C                     4 - subroutine will not change boundary faces
C                     8 - missing material interface and boundary
C                         faces created by the code will be removed
C                         from the final mesh
C                    16 - subroutine will not change surface points
C                   512 - tangled mesh will be fixed (testing)
C
C                    If only first two properties are required, set
C                    status = 5. Do not use the constant since status
C                    is the input/output parameter.
C
C
C     MaxSkipE  - the maximal number of skipped elements.
C                   The bad elements are put is a basket if their 
C                   cannot be fixed. These bad elements will be 
C                   skipped (temporary) from analysis. The basket 
C                   size is MaxSkipE. When it is full, all elements 
C                   are released back into the mesh.
C
C  ***BASIC 3: MaxQItr to iERR
C     MaxQItr   - the maximal number of local grid modifications
C
C     Metric(6,nP) - metric defined at mesh nodes. The metric is 
C                    a 3x3 symmetric matrix M_ij. Each column of
C                    this array stores the upper diagonal entries
C                    in the following order: M_11, M_22, M_33, 
C                    M_12, M_23, M_13.   
C
C     Quality   - target quality for the final grid
C                    (a positive number between 0 and 1)
C     rQuality  - actual quality of the final grid
C
C     MaxWr     - size of the working real*8 array
C     MaxWi     - size of the wirking integer array
C
C     rW(MaxWr) - real*8  working array. On output: 
C                   rW(1) - mesh generation time
C                   rW(2) - rQuality
C                   rW(3) - average size of elements with respect 
C                           to the given metrici (hStar in theory)
C                   rW(4:nE+3) - quality of elements
C
C
C     iW(MaxWi) - integer working array. On output:
C                 iW(1 : nP) - colors of mesh points as described
C                               in file color.fd
C                     T - tought point:   cannot be changed anyway
C                                         T-edges cannot be changed
C                     V - vertex:         cannot be changed anyway
C                    VB - edge point:     has to live on the edge
C                     B - surface point:  has to live on the surface
C                     I - internal point: free point inside the domain
C
C                   VBI - internal edge point:     similar to VB
C                    BI - internal boundary point: similar to  B
C
C                 iW(nP + 1) - the number of performed iterations
C
C     iPrint - nonnegative number 10 * A + B:
C                B - level of output information:
C                      0 - nothing 
C                      9 - maximum
C
C                A - output chanel number (between 0 and 99):
C                      A > 0 - the output goes to file aniMPI.log
C                      A = 0 - screen output is used
C
C
C     iERR   - error code:
C              0  - correct program termination
C           1000  - target quality was not reached (warning)
C           other - error as desribed in error.f 
C
C ================================================================
C
C  Note:
C       Input parameters:  MaxP, MaxF, MaxE, 
C                          nEStar, nPv, nFv, nEv, IPV, IFV, IEV, 
C                          flagAuto, MaxSkipE, MaxQItr,
C                          Metric, Quality, MaxWr, MaxWi, iPrint
C
C       Input / output:    nP, nF, nE, XYP, IPF, IPE,
C                          lbF, lbE, rW, iW
C
C       Output parameters: rQuality, iERR
C
C ================================================================
C
C  A possible choice for input parameters:
C       MaxP > nP
C       MaxF > nF
C       MaxE > nE
C       MaxSkipE = 200
C       MaxQItr = 25000
C       Quality = 0.3
C       MaxWr > 14 * MaxP + 2 * nP + MaxE  (approximate estimate)
C       MaxWi >  7 * MaxP + nP + 7 * MaxF + 18 * MaxE + 13 * nE
C                                                (approximately)
C
C ================================================================
C  A short version (mbaNodalShort) can be used when there are
C  no curvilinear boundaries. The following parameters are set
C  by default and have to be changed if necessary:
C
C       MaxSkipE =  200  
C       MaxQItr  =  min(50000, nE)  
C       Quality  =  0.3
C
C  Remark:
C    The reached mesh quality is returned in rW(2).
C
C ================================================================
C
C *** Authors: K. Lipnikov     (lipnikov@hotmail.com)
C              Yu. Vassilevski (vasilevs@dodo.inm.ras.ru)
C *** Date:    1997 - 2004
C *** Updates: see ChangeLog
C *** External routines: DSORT, DSYEV
C
C ==========================================================
C group (M)
C     Integer MaxP, MaxF, MaxE
      Real*8  XYP(3, *)
      Integer IPE(4, *), IPF(3, *)
      Integer lbF(*), lbE(*)

C group (Dev)
      Integer nPv, nFv, nEv, IPV(*), IFV(*), IEV(*) 
      Logical flagAuto
      Integer status

C group (Q)
C     Integer MaxSkipE, MaxQItr
      Real*8  Metric(6, *), Quality, rQuality
      Logical flagAnalytic

C group (W)
      Real*8  rW(*)
      Integer iW(*)

C group (Local variables)

      Logical flagFILE
      Integer chanelOUT, statusWork
      Real*8  hStar, tm0, tm1
      Character*80 message
      Character*(LineLenght) output(MaxLines)

      iERR = 0

      iPrintl   = mod(iPrint, 10)
      chanelOUT = iPrint / 10
      flagFILE  = iPrint.GE.10

      If(chanelOUT.GE.100) Call errMes(4201, 'ani_solution',
     &                          'output chanel number is wrong')


      If(flagFILE) Then
         Open(chanelOUT, file='aniMPI.log', status='UNKNOWN')
      End if


      nLines = 0
      nLoop = 1


c ... print header
      If(iPrintl.GE.1) Then
         If(flagFILE) Then
            Write(chanelOUT, 5001) Quality, nEStar, MaxQItr
         Else
            Write(*, 5001) Quality, nEStar, MaxQItr
         End if
      End if

      Call setStatus(flagAuto, status, iPrint)
      statusWork = status


c ... starting clocks
      Call mba3d_seconds(tm0)


c ... memory for sequantial running (is cleaned)
c     iW(1) is overloaded to return colors & nQItr
      iIPF = 1
      iIPE = iIPF + 4 * MaxF
      iIPP = iIPE + 5 * MaxE
      iICP = iIPP + MaxP
      iIHolP = iICP + MaxP
      iIHolF = iIHolP + MaxP
      iIHolE = iIHolF + MaxF
      iIEP = iIHolE + MaxE
      iIFE = iIEP + MaxP
      iIEE = iIFE + 4 * MaxE
      iL1E = iIEE + 4 * MaxE
      iL2E = iL1E + 2 * MaxE
      iIEPw = iL2E + MaxE
      inEPw = iIEPw + max(4 * MaxE, 3 * MaxF)
      iIPEw = inEPw + MaxP
      iiSE  = iIPEw + 4 * nE

      ilbP = 1

      miLINTRP = MaxWi - iiSE
      If(miLINTRP.LE.MaxP + MaxF) Then
         iERR = 1001
         Write(message,'(A,I10)')
     &        'The approximate size of iW is ', iiSE + MaxP + MaxF
         Call errMes(iERR, 'mbaNodal', message)
      End if


c ... memory for sequantial running (is cleaned)
c     rW(1) is overloaded to return total time
      iHesP = 1
      idG = iHesP + 6 * MaxP 
      iqE = idG + MaxP
      iXYPw = iqE + MaxE
      irSE = iXYPw + 3 * nP

      mrLINTRP = MaxWr - irSE
      If(mrLINTRP.LT.MaxE) Then
         iERR = 1002
         Write(message,'(A,I10)')
     &        'The approximate size of rW is ', irSE + MaxE
         Call errMes(iERR, 'mbaNodal', message)
      End if

      Do n = 1, MaxWr
         rW(n) = 0D0
      End do

      Do n = 1, MaxWi
         iW(n) = 0
      End do


c ... auxiliary structure of the data
      m = iIPF - 1
      Do n = 1, nF
         Do i = 1, 3
            iW(m + i) = IPF(i, n)
         End do
         iW(m + 4) = lbF(n)
         m = m + 4
      End do

      m = iIPE - 1
      Do n = 1, nE
         Do i = 1, 4
            iW(m + i) = IPE(i, n)
         End do
         iW(m + 5) = lbE(n)
         m = m + 5
      End do


c ... scale geometry to unit cube
      Call scale2Cube(nP, XYP, .TRUE.)


c ... setting the fixed metric for the future loops
      nPw = nP
      nEw = nE
      Call copyMeshData(nP, nE, XYP, Metric, IPE, 
     &                  rW(iXYPw), rW(iHesP), iW(iIPEw))

      Call checkMetric(nP, rW(iHesP))


c ... quality of tetrahedra is computed
      Call makQ(
     &     nLoop,
c group (M)
     &     nP, nE, XYP, IPE, nEv, IEV,
     &     nEStar, hStar,
c group (Q)
     &     rW(iHesP), rW(idG), rW(iqE))


c ... runing the basic algorithm for the global grid
      flagAnalytic = .FALSE.
      Call ani2(
c group (M)
     &     nP, MaxP, nF, MaxF, nE, MaxE, 
     &     XYP, iW(iIPF), iW(iIPE), 
     &     nEStar, hStar,
     &     iW(iICP), iW(iIPP), iW(iIEP),
     &     iW(iIFE), iW(iIEE),
     &     iW(iL1E), iW(iL2E),
     &     iW(iIHolP), iW(iIHolF), iW(iIHolE),
     &     iW(iIEPw), iW(inEPw),
     &     miLINTRP, mrLINTRP, iW(iIPEw), iW(iiSE),
     &     flagAuto, statusWork,
c group (Dev)
     &     nPv, nFv, nEv, IPV, IFV, IEV, 
c group (Q)
     &     MaxSkipE, MaxQItr, MaxBasketsGrid,
     &     rW(iHesP), Quality, rQuality, 
     &     rW(idG), rW(iqE), 
     &     nPw, nEw, rW(iXYPw), Metric, rW(irSE),
     &     MetricFunction_ani, flagAnalytic,
c group (ERR)
     &     flagFILE, chanelOUT, nLines, output,
     &     iPrintl, nQItr, iERR)

      If(iERR.NE.0 .AND. iERR.NE.1000)
     &   Call errMes(iERR, 'mbaNodal', 'memory problems')


      Call mba3d_seconds(tm1)
      tm1 = tm1 - tm0

      If(iPrintl.GE.1) Then
         If(flagFILE) Then
            Do i = 1, nLines
               Write(chanelOUT,'(A)') output(i)
            End do

            Write(chanelOUT, 5003) nQItr, rQuality, nP, nF, nE, tm1
         Else
            Write(*, 5003) nQItr, rQuality, nP, nF, nE, tm1
         End if
      End if


c ... close output file
      If(flagFILE) Close(chanelOUT)


c ... rescale geometry back
      Call scale2Cube(nP, XYP, .FALSE.)


c ... original structure of the data
      m = iIPF - 1
      Do n = 1, nF
         Do i = 1, 3
            IPF(i, n) = iW(m + i)
         End do
         lbF(n) = iW(m + 4)
         m = m + 4
      End do

      m = iIPE - 1
      Do n = 1, nE
         Do i = 1, 4
            IPE(i, n) = iW(m + i) 
         End do
         lbE(n) = iW(m + 5)
         m = m + 5
      End do

      Do n = 1, nP
         iW(n) = iW(iICP + n - 1)
      End do
      iW(nP + 1) = nQItr

      rW(1) = tm1
      rW(2) = rQuality
      rW(3) = hStar

      Do n = 1, nE
         rW(3 + n) = rW(iqE + n - 1)
      End do

      Return

 5001 Format(/,
     &    'STONE FLOWER! (1997-2008), version 2.1', /,
     &    'Target: Quality', F5.2, ' with', I9,
     &    ' tetrahedra for at most', I9, ' iterations',/)

 5002 Format(/,'Processor :', I4)

 5003 Format('Total:', I6, ' Q=', E10.4, '  #V,F,E:', I7,I8,I9,
     &       '  tm=', F6.1, 's',/)
      End Subroutine mbaNodal



C ================================================================
      Subroutine mbaNodalShort(
C ================================================================
c group (M)
     &      nP, MaxP, nF, MaxF, nE, MaxE,
     &      XYP, IPF, IPE, lbF,
     &      nEStar, status,
     &      Metric, 
c group (W)
     &      MaxWr, MaxWi, rW, iW,
     &      iPrint, iERR)
C ================================================================
      include 'status.fd'
      include 'magic.fd'
      include 'output.fd'
C ================================================================
C group (M)
      Real*8  XYP(3, *)
      Integer IPE(4, *), IPF(3, *)
      Integer lbF(*)
      Integer status

C group (Q)
      Real*8  Metric(6, *)

C group (W)
      Real*8  rW(*)
      Integer iW(*)

C group (Local variables)
      Real*8  Quality, rQuality
      Logical flagAuto
      Character*80 message

C ================================================================
      iERR = 0

      MaxWiAv =  7 * MaxP + nP + 9 * MaxF + 19 * MaxE + 13 * nE
      MaxWrAv = 14 * MaxP + 14 * nP + MaxE

      If(MaxWi.LE.MaxWiAv) Then
         iERR = 1001
         Write(message,'(A,I10)')
     &        'The approximate size of iW is ', MaxWiAv
         Call errMes(iERR, 'mbaNodalShort', message)
      End if

      If(MaxWr.LE.MaxWrAv) Then
         iERR = 1002
         Write(message,'(A,I10)')
     &        'The approximate size of rW is ', MaxWrAv
         Call errMes(iERR, 'mbaNodalShort', message)
      End if


c ... setup of missing parameter
      nPv = 0
      nFv = 0
      nEv = 0
      flagAuto = .TRUE.

      MaxSkipE = 200
      MaxQItr  = min(20 000, nE)
      MaxQItr  = max(1 000, MaxQItr)
      Quality  = 0.3D0


      iIPV = 1
      ilbE = iIPV + MaxF
      iIFV = ilbE + MaxE
      iIEV = iIFV + nFv + 1 
      iiEnd = iIEV + nEv + 1


      Do n = 0, nE - 1
         iW(ilbE + n) = 1
      End do

      Call mbaNodal(
c group (M)
     &      nP, MaxP, nF, MaxF, nE, MaxE,
     &      XYP, IPF, IPE, lbF, iW(ilbE),
     &      nEStar, 
c group (Dev)
     &      nPv, nFv, nEv, iW(iIPV), iW(iIFV), iW(iIEV), 
     &      flagAuto, status,
c group (Q)
     &      MaxSkipE, MaxQItr,
     &      Metric, Quality, rQuality, 
c group (W)
     &      MaxWr, MaxWi - iiEnd, rW, iW(iiEnd),
     &      iPrint, iERR)

      Return
      End Subroutine mbaNodalShort
C
      End Module mba3d_mba_nodal
