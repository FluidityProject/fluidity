      module mba2d_module

      contains

C ==========================================================
      Subroutine mbaNodal(
C ==========================================================
c group (M)
     &      nP, MaxP, nF, MaxF, nE, MaxE, nPv,
     &      XYP, IPF, IPE, IPV,
     &      CrvFunction, ParCrv, iFnc,
     &      nEStar, 
c group (D)
     &      nFv, nEv, IFV, IEV, lbE,
     &      flagAuto, status,
c group (Q)
     &      MaxSkipE, MaxQItr,
     &      Metric, Quality, rQuality, 
c group (W)
     &      MaxWr, MaxWi, rW, iW,
     &      iPrint, iERR)
C ==========================================================
      include 'lintrp.fd'
      include 'status.fd'
C ==========================================================
C This routine adapts the mesh to the discrete metric given 
C at mesh nodes.
C ==========================================================
C  VARIABLES & PARAMETER:
C  
C     nP, MaxP  - the real and maximal number of points (P)
C     nF, MaxF  - the real and maximal number of boundary and inner
C     edges (F)
C     nE, MaxE  - the real and maximal number of triangles (E)
C     nPv       - the number of fixed points (V)
C
C     XYP(2, MaxP) - the Cartesian coordinates of mesh points
C     IPE(3, MaxE) - connectivity list of triangles
C     IPF(4, MaxF) - column 1 &  2 - connectivity list of boundary and
C     inner edges
C                    column 3      - number in the list of
C                    parametrizations:
C
C                                   =0     - this edge is the linear
C                                   segment
C                                   =n > 0 - ParCrv(*, n) gives a
C                                   parametrization
C                                            of this edge and iFnc(n)
C                                            gives
C                                            a function number for
C                                            computing the 
C                                            Cartesian coordinates (see
C                                            calCrv())
C                    column 4      - boundary identificator 
C                                    (example: unit square has 4
C                                    boundaries which
C                                     may have or not different
C                                     identificators.
C                                     In order to automatically
C                                     recognize corners
C                                     of the square, boundaries have to
C                                     have 
C                                     different colors. It is not
C                                     required if 
C                                     the corner points are in the list
C                                     of fixed
C                                     points. See colors.fd for more
C                                     details.)
C
C     IPV(nPv)     - list of fixed points
C
C
C     CrvFunction  - user-created routine that computes the 
C                    Cartesian coordinates of a point xyc from 
C                    its parametric coordinate tc:
C
C                    Subroutine CrvFunction(tc, xyc, iFnc)
C
C                      tc     - parametric coordinate of point xyc
C                      xyc(2) - Cartesian coordinate of the same point
C                      iFnc   - the function number associated with a
C                               curved edge
C
C                      On input :  tc, iFnc
C                      On output:  xyc(2)
C
C     ParCrv(2, MaxF) - linear parametrizations of curvilinear edges 
C                       column 1 - parameter for the starting point
C                       column 2 - parameter for the terminal point
C
C                       parameters for the inner points are computed by
C                       the
C                       linear interpolation between two given numbers
C
C                       the Cartesian coordinates are computed by
C                       user-given formulas defined inside CrvFunction.
C
C     iFnc(MaxF)      - function number for computing the Cartesian
C     coordinates
C
C     nEstar - the desired number of triangles
C
C     
C     nFv       - number of fixed edges
C     nEv       - number of fixed triangles
C     IFV(nFv)  - list of fixed edges
C     IEV(nEv)  - list of fixed triangles
C
C     lbE(MaxE) - element indentificator (a positive number)
C
C     flagAuto  - flag controling the mesh generation:
C                 TRUE  - automatic recovering of missing mesh elements
C                 FALSE - rigorous checking of input data 
C
C     status    - advanced control of mesh generation:
C                 0 or negative - no additional control
C                 positive      - enforce some of the mesh properties
C                      The detailed description of available properties
C                      is in file status.fd. Variable status is equal 
C                      to the sum of positive numbers corresponding
C                      to the desired properties. Here is the list of 
C                      user-controled properties:
C  
C                      1  - the final mesh will not have boundary
C                      elements;
C                      4  - the algorithm will not change boundary
C                      edges;
C                      8  - the missing material interfaces and
C                      boundary;
C                           edges created by the code will be removed 
C                           from the final mesh;
C                      16 - the algorithm will not change boundary
C                      points.
C
C                      If only two first properties are required, set 
C                      status = 5. Do not use tne number since status
C                      is input/output parameter.
C
C
C     MaxSkipE  - the maximal number of skipped triangles
C     MaxQItr   - the maximal number of local grid modifications
C
C     Metric(3, nP) - real array containing the metric defined at
C                     mesh points. The metric is a 2x2 positive 
C                     definite symmetric tensor:
C
C                            M11   M12   
C                   Metric = 
C                            M12   M22
C
C                   Each column of this array stores the upper
C                   triangular entries in the following order:
C                   M11, M22, and M12.
C
C     Quality   - the requested quality for the final grid
C                 (a positive number between 0 and 1)
C     rQuality  - the reached quality for the final grid
C
C     MaxWr     - the maximal space for real arrays
C     MaxWi     - the maximal space for integer arrays
C
C     rW(MaxWr)  - the real  working array. On output: 
C                  rW(1) - mesh generation time
C                  rW(2) - rQuality
C                  rW(3) - average size of triangles (with respect to 
C                          the given metric) [hStar]
C
C     iW(MaxWi)  - the integer working array. On output:
C                  iW(1 : nP) - colors of mesh points as 
C                               described in file color.fd
C                  iW(nP + 1) - the number of performed iterations
C       
C                  If memory allocation unsufficient, the output
C                  is different:
C                  iW(1) - the required integer memory allocation
C                  iW(2) - the required real  memory allocation
C
C
C     iPrint - nonnegative number 10 * A + B:
C
C              B  - the level of output information:
C                   0 - nothing
C                   9 - maximum
C
C              A - the output chanel number (between 0 and 99):
C                  A > 0 - the output goes to file aniMPI.log
C                  A = 0 - screen output is used
c
C     iERR - the error code:
C             0 - a correct termination
C          1000 - the quality has not been reached
C
C ================================================================
C
C  Note:
C       Input parameters:  MaxP, MaxF, MaxE, nPv,
C                          IPV, IFV, IEV, lbE, flagAuto, 
C                          nEStar, MaxSkipE, MaxQItr,
C                          Quality, MaxWr, MaxWi, iPrint
C
C       Input / output:    nP, nF, nE, XYP, IPF, IPE,
C                          ParCrv, iFnc, Sol, status, rW, iW
C
C       Output parameters: rQuality, iERR
C
C ================================================================
C
C  A possible choice for input parameters:
C       MaxP > nP
C       MaxF > nF
C       MaxE > nE
C       nPv >= 0
C       nFv >= 0
C       nEv >= 0
C       MaxSkipE = 100
C       MaxQItr = 15000
C       Quality = 0.6
C       MaxWr > 4 * MaxP + 10 * nP + MaxF + MaxE (approximate estimate)
C       MaxWi > 6 * MaxP + 10 * nP + 19 * MaxF + 11 * MaxE + 12 * nE
C       (approximately)
C
C *** Authors: K. Lipnikov (lipnikov@hotmail.com)
C              Y. Vassilevski (vasilevs@dodo.inm.ras.ru)
C ==========================================================
C group (M)
      Integer MaxP, MaxF, MaxE
      real  XYP(2, *)
      Integer IPE(3, *), IPF(4, *), IPV(*)

      EXTERNAL CrvFunction
      real   ParCrv(2, *)
      Integer  iFnc(*)

C group (D)
      Integer nFv, nEv
      Integer IFV(*), IEV(*), lbE(*)
      
      Logical flagAuto
      Integer status

C group (Q)
      Integer MaxSkipE, MaxQItr
      real  Metric(3, *)
      real  Quality, rQuality
      Logical flagAnalytic

C group (W)
      real  rW(*)
      Integer iW(*)

C group (Local variables)
      real  hStar

      Integer  MetricFunction_ani
      EXTERNAL MetricFunction_ani

C ==========================================================
C group (Common blocks)
      Integer iDomBnd, iMatBnd
      Common /aniBND/ iDomBnd, iMatBnd
 
      real  refXYP(2), scaXYP(2)
      Common /rescale/refXYP, scaXYP

C ==========================================================
      iERR = 0

      Do i = 1, 2
         refXYP(i) = 0D0
         scaXYP(i) = 1D0
      End do


c ... refine initial mesh when nE is very small
c ... it increases robustness of the code 
      Do while(nE < nEStar / 15 .AND. nE.LE.500 .AND. nEv+nFv.EQ.0)
         iIFE = 1
         iiW  = iIFE + 3 * nE
         nWi  = iiW  + 3 * nE + nP 
         If(nWi.GT.MaxWi) goto 100

         If(iPrint.GE.1) Write(*,5001) nP, nE

         Call uniformRefinement(
     &        nP, MaxP, nF, MaxF, nE, MaxE,
     &        XYP, IPF, IPE, lbE,
     &        CrvFunction, ParCrv, iFnc, iW(iIFE),
     &        Metric, 3, iW(iiW), MaxWi)
      End do 


 100  miLINTRP = 10 * nP + 3 * nE + 6
      mrLINTRP =  4 * nP + MaxH + 4

      inEt = 1
      inStept = inEt + MaxF
      inL2t = inStept + 4 * MaxF
      iLFnc = inL2t + MaxF
      iILt  = iLFnc + MaxF
      iL1Et = iILt + MaxF
      iL2Et = iL1Et + 2 * MaxF
      iIHolP = iL2Et + 2 * MaxF
      iIHolF = iIHolP + MaxP
      iIHolE = iIHolF + MaxF
      iICP = iIHolE + MaxE
      iIEP = iICP + MaxP
      iIFE = iIEP + MaxP
      iIEE = iIFE + 3 * MaxE
      iL1E = iIEE + 3 * MaxE
      iL2E = iL1E + 2 * MaxE
      iIPEw = iL2E + 2 * MaxE
      iiSE  = iIPEw + 3 * nE
      iIEPw = iiSE + miLINTRP
c ... we need twice less memory for backReferences
      inEPw = iIEPw + max(6 * nE, 4 * MaxF)
      nWi   = inEPw + max(3 * MaxP, 2 * MaxF)


      itE = 1
      idG = itE + MaxF
      iqE = idG + MaxP
      iHesPw = iqE + MaxE
      iXYPw = iHesPw + 3 * nP
      irSE = iXYPw + 2 * nP
      nWr  = irSE + max(mrLINTRP, max(nE, MaxF))


      iW(1) = nWi
      iW(2) = nWr
      If(nWi.GT.MaxWi) Then
         iERR = 1001
         goto 1000
      End if

      If(nWr.GT.MaxWr) Then
         iERR = 1002
         goto 1000
      End if


      Do n = 1, nWr
         rW(n) = 0D0
      End do

      Do n = 1, nWi
         iW(n) = 0
      End do


c ... scale geometry to unit cube
      Call scale2Square(nP, XYP, .TRUE.)


c ... print Ani2D header
      If(iPrint.GE.1) Write(*, 5004) Quality, nEStar, MaxQItr


c ... set up default status
      Call setStatus(flagAuto, status, iPrint)


c ... call the main module
      flagAnalytic = .FALSE.
      Call ani2(
c group (M)
     &      nP, MaxP, nF, MaxF, nE, MaxE, nPv,
     &      XYP, IPF, IPE, IPV,
     &      CrvFunction, ParCrv, iFnc,
     &      nEStar, hStar,
     &      iW(iICP), iW(iIEP),
     &      iW(iIFE), iW(iIEE),
     &      iW(iL1E), iW(iL2E),
     &      iW(iIHolP), iW(iIHolF), iW(iIHolE),
     &      iW(iIEPw), iW(inEPw),
     &      iW(iIPEw), iW(iiSE),
c group (Dev)
     &      nFv, nEv, IFV, IEV, lbE,
     &      flagAuto, status,
c group (CRV)
     &      iW(iL1Et), iW(iL2Et), rW(itE),
     &      iW(inL2t), iW(inStept), iW(inEt),
     &      iW(iLFnc), iW(iILt),
c group (Q)
     &      MaxSkipE, MaxQItr, nQItr,
     &      Metric, Quality, rQuality,
     &      rW(idG), rW(iqE), rW(iXYPw), rW(iHesPw), rW(irSE),
     &      MetricFunction_ani, flagAnalytic,
c group (ERR)
     &      iPrint, iERR)


c ... rescale geometry back
      Call scale2Square(nP, XYP, .FALSE.)


c ... returning sadditional information
      Do n = 1, nP
         iW(n) = iW(iICP + n - 1)
      End do
      iW(nP + 1) = nQItr

      rW(1) = 0D0
      rW(2) = rQuality
      rW(3) = hStar


 1000 If(iERR.EQ.0 .OR. iERR.EQ.1000) Return
      Call errMes(iERR, 'mbaMetric', 
     &            'See error.f for error description')

      Return

 5001 Format('Refining mesh:', I6, ' pts  and', I7, ' elements')

 5004 Format(/,
     &    'STONE FLOWER! (1997-2007), version 2.0', /,
     &    'Target: Quality', F5.2, ' with', I8, 
     &    ' triangles for at most', I8, ' iterations') 
      End subroutine


      end module mba2d_module
