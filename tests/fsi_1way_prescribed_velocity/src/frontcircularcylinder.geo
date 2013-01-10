// lc= 0.015;
lc= 0.0025;
mc= 0.01;

// Create circle in middle of cylinder:
Point(1) = {0.2,0.2,0.0,mc};
Point(2) = {0.25,0.2,0.0,lc};
Point(3) = {0.2,0.25,0.0,lc};
Circle(1) = {2,1,3};
Line(101) = {1,2};

Point(4) = {0.15,0.2,0.0,lc};
Point(5) = {0.2,0.15,0.0,lc};
Circle(2) = {3,1,4};
Line(102) = {1,3};
Circle(3) = {4,1,5};
Line(103) = {1,4};
Circle(4) = {5,1,2};
Line(104) = {1,5};
//Create surface of middle circle of cylinder:
//Line Loop(5) = {4, 1, 2, 3};
Line Loop(5) = {101, 1, -102};
Line Loop(6) = {102, 2, -103};
Line Loop(7) = {103, 3, -104};
Line Loop(8) = {104, 4, -101};
Ruled Surface(11) = {5};
Ruled Surface(12) = {6};
Ruled Surface(13) = {7};
Ruled Surface(14) = {8};

//Physical lines:
Physical Line(1) = {2, 1, 4, 3};
//Physical surface:
Physical Surface(2) = {11, 12, 13, 14};


