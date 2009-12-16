!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    C.Pain@Imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA
module RGPTWE_module
  implicit none
contains

  REAL FUNCTION RGPTWE(IG,ND,WEIGHT)
    IMPLICIT NONE
    !     NB If WEIGHT is TRUE in function RGPTWE then return the Gauss-pt weight 
    !     else return the Gauss-pt. 
    !     NB there are ND Gauss points we are looking for either the 
    !     weight or the x-coord of the IG'th Gauss point. 
    INTEGER IG,ND
    LOGICAL WEIGHT

    IF(WEIGHT) THEN
       GO TO (10,20,30,40,50,60,70,80,90,100) ND
       !     +++++++++++++++++++++++++++++++
       !     For N=1 +++++++++++++++++++++++
10     CONTINUE
       RGPTWE=2.0
       GO TO 1000
       !     For N=2 +++++++++++++++++++++++
20     CONTINUE
       RGPTWE=1.0
       GO TO 1000
       ! For N=3 +++++++++++++++++++++++
30     CONTINUE
       GO TO (11,12,11) IG
11     RGPTWE= 0.555555555555556
       GO TO 1000
12     RGPTWE= 0.888888888888889
       GO TO 1000
       ! For N=4 +++++++++++++++++++++++
40     CONTINUE
       GO TO (21,22,22,21) IG
21     RGPTWE= 0.347854845137454
       GO TO 1000
22     RGPTWE= 0.652145154862546
       GO TO 1000
       ! For N=5 +++++++++++++++++++++++
50     CONTINUE
       GO TO (31,32,33,32,31) IG
31     RGPTWE= 0.236926885056189
       GO TO 1000
32     RGPTWE= 0.478628670499366
       GO TO 1000
33     RGPTWE= 0.568888888888889
       GO TO 1000
       ! For N=6 +++++++++++++++++++++++
60     CONTINUE
       GO TO (41,42,43,43,42,41) IG
41     RGPTWE= 0.171324492379170
       GO TO 1000
42     RGPTWE= 0.360761573048139
       GO TO 1000
43     RGPTWE= 0.467913934572691
       GO TO 1000
       ! For N=7 +++++++++++++++++++++++
70     CONTINUE
       GO TO (51,52,53,54,53,52,51) IG
51     RGPTWE= 0.129484966168870
       GO TO 1000
52     RGPTWE= 0.279705391489277
       GO TO 1000
53     RGPTWE= 0.381830050505119
       GO TO 1000
54     RGPTWE= 0.417959183673469
       GO TO 1000
       ! For N=8 +++++++++++++++++++++++
80     CONTINUE
       GO TO (61,62,63,64,64,63,62,61) IG
61     RGPTWE= 0.101228536290376
       GO TO 1000
62     RGPTWE= 0.222381034453374
       GO TO 1000
63     RGPTWE= 0.313706645877877
       GO TO 1000
64     RGPTWE= 0.362683783378362
       GO TO 1000
       ! For N=9 +++++++++++++++++++++++
90     CONTINUE
       GO TO (71,72,73,74,75,74,73,72,71) IG
71     RGPTWE= 0.081274388361574
       GO TO 1000
72     RGPTWE= 0.180648160694857
       GO TO 1000
73     RGPTWE= 0.260610696402935
       GO TO 1000
74     RGPTWE= 0.312347077040003
       GO TO 1000
75     RGPTWE= 0.330239355001260
       GO TO 1000
       ! For N=10 +++++++++++++++++++++++
100    CONTINUE
       GO TO (81,82,83,84,85,85,84,83,82,81) IG
81     RGPTWE= 0.066671344308688
       GO TO 1000
82     RGPTWE= 0.149451349150581
       GO TO 1000
83     RGPTWE= 0.219086362515982
       GO TO 1000
84     RGPTWE= 0.269266719309996
       GO TO 1000
85     RGPTWE= 0.295524224714753
       !
1000   CONTINUE
    ELSE
       GO TO (210,220,230,240,250,260,270,280,290,200) ND
       ! +++++++++++++++++++++++++++++++
       ! For N=1 +++++++++++++++++++++++ THE GAUSS POINTS...
210    CONTINUE
       RGPTWE=0.0
       GO TO 2000
       ! For N=2 +++++++++++++++++++++++
220    CONTINUE
       RGPTWE= 0.577350269189626
       GO TO 2000
       ! For N=3 +++++++++++++++++++++++
230    CONTINUE
       GO TO (211,212,211) IG
211    RGPTWE= 0.774596669241483
       GO TO 2000
212    RGPTWE= 0.0
       GO TO 2000
       ! For N=4 +++++++++++++++++++++++
240    CONTINUE
       GO TO (221,222,222,221) IG
221    RGPTWE= 0.861136311594953
       GO TO 2000
222    RGPTWE= 0.339981043584856
       GO TO 2000
       ! For N=5 +++++++++++++++++++++++
250    CONTINUE
       GO TO (231,232,233,232,231) IG
231    RGPTWE= 0.906179845938664
       GO TO 2000
232    RGPTWE= 0.538469310105683
       GO TO 2000
233    RGPTWE= 0.0
       GO TO 2000
       ! For N=6 +++++++++++++++++++++++
260    CONTINUE
       GO TO (241,242,243,243,242,241) IG
241    RGPTWE= 0.932469514203152
       GO TO 2000
242    RGPTWE= 0.661209386466265
       GO TO 2000
243    RGPTWE= 0.238619186083197
       GO TO 2000
       ! For N=7 +++++++++++++++++++++++
270    CONTINUE
       GO TO (251,252,253,254,253,252,251) IG
251    RGPTWE= 0.949107912342759
       GO TO 2000
252    RGPTWE= 0.741531185599394
       GO TO 2000
253    RGPTWE= 0.405845151377397
       GO TO 2000
254    RGPTWE= 0.0
       GO TO 2000
       ! For N=8 +++++++++++++++++++++++
280    CONTINUE
       GO TO (261,262,263,264,264,263,262,261) IG
261    RGPTWE= 0.960289856497536
       GO TO 2000
262    RGPTWE= 0.796666477413627
       GO TO 2000
263    RGPTWE= 0.525532409916329
       GO TO 2000
264    RGPTWE= 0.183434642495650
       GO TO 2000
       ! For N=9 +++++++++++++++++++++++
290    CONTINUE
       GO TO (271,272,273,274,275,274,273,272,271) IG
271    RGPTWE= 0.968160239507626
       GO TO 2000
272    RGPTWE= 0.836031107326636
       GO TO 2000
273    RGPTWE= 0.613371432700590
       GO TO 2000
274    RGPTWE= 0.324253423403809
       GO TO 2000
275    RGPTWE= 0.0
       GO TO 2000
       ! For N=10 +++++++++++++++++++++++
200    CONTINUE
       GO TO (281,282,283,284,285,285,284,283,282,281) IG
281    RGPTWE= 0.973906528517172
       GO TO 2000
282    RGPTWE= 0.865063366688985
       GO TO 2000
283    RGPTWE= 0.679409568299024
       GO TO 2000
284    RGPTWE= 0.433395394129247
       GO TO 2000
285    RGPTWE= 0.148874338981631
       !
2000   CONTINUE
       IF(IG.LE.INT((ND/2)+0.1)) RGPTWE=-RGPTWE
    ENDIF
  END FUNCTION RGPTWE
end module RGPTWE_module
