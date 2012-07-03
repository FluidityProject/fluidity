Point(1) = {0.0, 0.0, 0.0, 40};
Point(2) = {0.0, 100.0, 0.0, 40};
Point(3) = {0.0, 1000.0, 0.0, 40};
Point(4) = {5000.0, 1000.0, 0.0, 40};
Point(5) = {5000.0, 0.0, 0.0, 40};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,1};
Line Loop(6) = {1,2,3,4,5};

Plane Surface(7) = {6};

// Top of the box
Physical Line(333) = {3};
// Box bottom
Physical Line(444) = {5};
// Box left side - the rest of it
Physical Line(555) = {2};
// Box right side
Physical Line(666) = {4};
// Box left side - inflow
Physical Line(999) = {1};

// This is just to ensure all the interior
// elements get written out. 
Physical Surface(10) = {7};
