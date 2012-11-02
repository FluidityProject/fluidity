H = 0.1;
x1 = 3*H;
x2 = 6*H;
h = 2*H;
b = 7*H;
dx = H/3;

// Points of outer (upper half-)cube
Point(1) = {-x1, 0, b/2, dx};
Point(2) = {-x1, h, b/2, dx};
Point(3) = {H+x2, 0, b/2, dx};
Point(4) = {H+x2, h, b/2, dx};
Point(5) = {-x1, 0, 0, dx};
Point(6) = {-x1, h, 0, dx};
Point(7) = {H+x2, 0, 0, dx};
Point(8) = {H+x2, h, 0, dx};

// Points of inner (upper half-)cube
Point(9) = {0, 0, H/2, dx};
Point(10) = {0, H, H/2, dx};
Point(11) = {H, 0, H/2, dx};
Point(12) = {H, H, H/2, dx};
Point(13) = {0, 0, 0, dx};
Point(14) = {0, H, 0, dx};
Point(15) = {H, 0, 0, dx};
Point(16) = {H, H, 0, dx};

Line(1) = {1,2};
Line(2) = {1,5};
Line(3) = {5,6};
Line(4) = {6,2};
Line(5) = {2,4};
Line(6) = {1,3};
Line(7) = {4,3};
Line(8) = {4,8};
Line(9) = {8,6};
Line(10) = {7,3};
Line(11) = {7,15};
Line(12) = {15,16};
Line(13) = {16,14};
Line(14) = {14,13};
Line(15) = {13,5};
Line(16) = {7,8};

Line(17) = {16,12};
Line(18) = {12,10};
Line(19) = {10,9};
Line(20) = {11,9};
Line(21) = {15,11};
Line(22) = {13,9};
Line(23) = {12,11};
Line(24) = {10,14};

Line Loop(1) = {1,-4,-3,-2}; // Inlet
Line Loop(2) = {6,-7,-5,-1}; // Side along z = b/2
Line Loop(3) = {3,-9,-16,11,12,13,14,15}; // Side along z = 0
Line Loop(4) = {7,-10,16,-8}; // Outlet
Line Loop(5) = {5,8,9,4}; // Top
Line Loop(6) = {-2,6,-10,11,21,20,-22,15}; // Bottom
Line Loop(7) = {22,-19,24,14}; // Inner cube sides along x = 0
Line Loop(8) = {21,-23,-17,-12}; // Inner cube sides along x = H
Line Loop(9) = {24,-13,17,18}; // Inner cube top
Line Loop(10) = {20,-19,-18,23}; // Inner cube sides along z = H/2

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};
Plane Surface(7) = {7};
Plane Surface(8) = {8};
Plane Surface(9) = {9};
Plane Surface(10) = {10};

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
Physical Line(5) = {5};
Physical Line(6) = {6};
Physical Line(7) = {7};
Physical Line(8) = {8};
Physical Line(9) = {9};
Physical Line(10) = {10};

// This is just to ensure all the interior
// elements get written out. 
//Physical Surface(11) = {1,2,3,4,5,6,7,8,9,10};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
Physical Surface(5) = {5};
Physical Surface(6) = {6};
Physical Surface(7) = {7};
Physical Surface(8) = {8};
Physical Surface(9) = {9};
Physical Surface(10) = {10};

Surface Loop(800) = {1,2,3,4,5,6,7,8,9,10};
Volume(900) = {800};
Physical Volume(999) = {900};