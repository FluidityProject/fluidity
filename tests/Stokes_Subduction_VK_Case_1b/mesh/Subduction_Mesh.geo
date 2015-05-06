// Gmsh project created on Fri Dec  4 14:59:47 2009
// We reproduce a subduction zone on a (660x600) km domain. 
// The kinematically moving slab (see image, on the left) subducts under a 50km-thick plate (up, right) and hence dynamically inducing a flow in the upper-mantle wedge (down, on the right). 
// We therefore are trying to match the benchmark results established by vanKeken2008 for this problem, computing the thermal field throughout the domain by use of finite element methods. 
// The structured, non-adaptive mesh shown in the images contains for instance the sought information at the mesh nodes, spaced out 2.5 km from each other near the critical point of contact between the three// different regions of our domain, and 7.5km everywhere else. In this particular case, the test uses a prescribed wedge flow (Batchelor1967) and is mainly intended to check for the validity of initial condi// tions at the various boundaries. 


Point(1) = {0,0,0,2000};
Point(2) = {660000,0,0,10000};
Point(3) = {660000,-50000,0,10000};
Point(4) = {660000,-600000,0,10000};
Point(5) = {600000,-600000,0,10000};
Point(6) = {0,-600000,0,10000};
Point(7) = {50000,-50000,0,2000};
Point(8) = {660000,-420000,0,10000};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,8};
Line(4) = {4,5};
Line(5) = {5,6};
Line(7) = {1,7};
Line(9) = {7,3};
Line(10) = {8,4};
Line(11) = {7,5};
Line(14) = {6,1};

Physical Line(10) = {1};
Physical Line(11) = {2};
Physical Line(12) = {3};
Physical Line(13) = {4};
Physical Line(19) = {5};
Physical Line(15) = {14};
Physical Line(16) = {10};
// Physical Line(17) = {9};
// Physical Line(18) = {11};
// Physical Line(20) = {7};

Line Loop(22) = {1,2,-9,-7};
Plane Surface(23) = {22};
Physical Surface(24) = {23};
Line Loop(25) = {9,3,10,4,-11};
Plane Surface(26) = {25};
Physical Surface(27) = {26};
Line Loop(28) = {5,14,7,11};
Plane Surface(29) = {28};
Physical Surface(30) = {29};

Physical Line(31) = {11};
Physical Line(32) = {9};
