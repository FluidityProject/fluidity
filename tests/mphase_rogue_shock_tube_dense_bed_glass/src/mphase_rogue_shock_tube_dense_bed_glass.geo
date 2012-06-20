Point(1) = {0.0, 0.0, 0.0, 0.005};
Point(2) = {0.0, 7.0, 0.0, 0.005};
Point(3) = {0.13, 7.0, 0.0, 0.005};
Point(4) = {0.13, 0.0, 0.0, 0.005};
Line(1) = {3,4};
Line(2) = {4,1};
Line(3) = {1,2};
Line(4) = {2,3};
Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};
// Top of the box
Physical Line(333) = {4};
// Box sides
Physical Line(666) = {3,1};
// Box bottom
Physical Line(444) = {2};
// This is just to ensure all the interior
// elements get written out. 
Physical Surface(10) = {6};
