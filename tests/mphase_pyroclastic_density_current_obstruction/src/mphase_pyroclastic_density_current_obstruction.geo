Point(1) = {0.0, 0.0, 0.0, 40};
Point(2) = {0.0, 100.0, 0.0, 40};
Point(3) = {0.0, 1000.0, 0.0, 40};
Point(4) = {5000.0, 1000.0, 0.0, 40};
Point(5) = {5000.0, 0.0, 0.0, 40};
Point(6) = {750.0, 0.0, 0.0, 40};
Point(7) = {750.0, 75.0, 0.0, 40};
Point(8) = {500.0, 75.0, 0.0, 40};
Point(9) = {500.0, 0.0, 0.0, 40};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,9};
Line(9) = {9,1};
Line Loop(10) = {1,2,3,4,5,6,7,8,9};

Plane Surface(11) = {10};

// Top of the box
Physical Line(333) = {3};
// Box bottom
Physical Line(444) = {5,9};
// Box left side - the rest of it
Physical Line(555) = {2};
// Box right side
Physical Line(666) = {4};
// Box left side - inflow
Physical Line(999) = {1};
// Obstruction - vertical sides
Physical Line(111) = {6,8};
// Obstruction - top
Physical Line(222) = {7};

// This is just to ensure all the interior
// elements get written out. 
Physical Surface(12) = {11};
