lc =0.01;
Point(1) = {0.0, 0.0, 0.0, lc};
Point(2) = {0.0, 0.45, 0.0, lc};
Point(3) = {0.2, 0.45, 0.0, lc};
Point(4) = {0.2, 0.0, 0.0, lc};
Point(5) = {0.11, 0.0, 0.0, lc};
Point(6) = {0.09, 0.0, 0.0, lc};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};
Line Loop(5) = {1,2,3,4,5,6};
Plane Surface(6) = {5};
// Bottom of the box part 1
Physical Line(331) = {6};
// Bottom of the box part 2
Physical Line(332) = {5};
// Bottom of the box part 3
Physical Line(333) = {4};
// Box sides left
Physical Line(661) = {1};
// Box sides right
Physical Line(662) = {3};
// Box top
Physical Line(444) = {2};
// This is just to ensure all the interior
// elements get written out. 
Physical Surface(10) = {6};
