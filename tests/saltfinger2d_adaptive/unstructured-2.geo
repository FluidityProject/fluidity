Point(1) = {0,0,0,0.004};
Point(2) = {0,2,0,0.004};
Point(3) = {4,2,0,0.004};
Point(4) = {4,0,0,0.004};Line(1) = {2,1};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {2,3,4,-1};
Plane Surface(6) = {5};
//top
Physical Line(7) = {2};
//bottom
Physical Line(8) = {4};
//left
Physical Line(9) = {1};
//right
Physical Line(10) = {3};

Physical Surface(11) = {6};
