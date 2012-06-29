Point(1) = {0.0, 0.0, 0.0, 1};
Point(2) = {75, 0.0, 0.0, 1};
Point(3) = {80, 0.0, 0.0, 1};
Point(4) = {140.0, 0.0, 0.0, 1};
Point(5) = {140.0, 50.0, 0.0, 1};
Point(6) = {0.0, 50.0, 0.0, 1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};
Line Loop(7) = {1,2,3,4,5,6};

Plane Surface(8) = {7};

// Box bottom - vent
Physical Line(111) = {2};
// Box bottom - the rest of it
Physical Line(222) = {1,3};
// Top of the box
Physical Line(333) = {5};

// Box left side
Physical Line(888) = {6};
// Box right side
Physical Line(999) = {4};

// This is just to ensure all the interior
// elements get written out. 
Physical Surface(10) = {8};
