// element size
s=0.075; // coarse
//s=0.025; // fine

Point(1) = {0. , 0. , 0., s};
Point(2) = {0.5, 0. , 0., s};
Point(3) = {1. , 0. , 0., s};
Point(4) = {0. , 0.5, 0., s};
Point(5) = {0.5, 0.5, 0., s};
Point(6) = {1. , 0.5, 0., s};
Point(7) = {0. , 1. , 0., s};
Point(8) = {0.5, 1. , 0., s};
Point(9) = {1. , 1. , 0., s};


Point(11) = {0.25 , 0.25 , 0., s};
Point(12) = {0.5 , 0.25 , 0., s};
Point(13) = {0.75 , 0.25 , 0., s};
Point(14) = {0.25 , 0.5 , 0., s};
Point(15) = {0.5 , 0.5 , 0., s};
Point(16) = {0.75 , 0.5 , 0., s};
Point(17) = {0.25 , 0.75 , 0., s};
Point(18) = {0.5 , 0.75 , 0., s};
Point(19) = {0.75 , 0.75 , 0., s};

Line(1) = {1, 3};
Line(2) = {3, 13};
Line(3) = {13, 11};
Line(4) = {11, 1};


Line(5) = {1, 7};
Line(6) = {7, 17};
Line(7) = {17, 11};

Line(8) = {17, 19};
Line(9) = {19, 9};
Line(10) = {9, 7};

Line(11) = {13, 19};
Line(12) = {9, 3};


Line Loop(111) = {1, 2,  3, 4};
Plane Surface(112) = {111};

Line Loop(113) = {4, 5,  6, 7};
Plane Surface(114) = {113};

Line Loop(115) = {9, 10,  6, 8};
Plane Surface(116) = {115};

Line Loop(117) = {2, 11,  9, 12};
Plane Surface(118) = {117};

Line Loop(119) = {-3, 11,  -8, 7};
Plane Surface(120) = {119};


// inflow
Physical Line(21) = {5};
// outflow
Physical Line(22) = {12};
// sides
Physical Line(23) = {1, 10};
// surface
Physical Surface(24) = {112, 114, 116, 118};
Physical Surface(25) = {120};
