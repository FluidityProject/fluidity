Point(1) = {0.0, 0.0, 0.0, 0.03};
Point(2) = {0.3, 0.0, 0.0, 0.03};
Point(3) = {0.3, 0.0, 0.3, 0.03};
Point(4) = {0.0, 0.0, 0.3, 0.03};
Point(5) = {0.0, 0.7, 0.0, 0.03};
Point(6) = {0.3, 0.7, 0.0, 0.03};
Point(7) = {0.3, 0.7, 0.3, 0.03};
Point(8) = {0.0, 0.7, 0.3, 0.03};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,5};

Line(9) = {5,1};
Line(10) = {6,2};
Line(11) = {7,3};
Line(12) = {8,4};

// Bottom surface
Line Loop(5) = {1,2,3,4};

// Side surface, nearest to the x axis
Line Loop(6) = {1,-10,-5,9};

// Side surface, furthest from the x axis
Line Loop(7) = {-3,-11,7,12};

// Side surface, nearest to the z axis
Line Loop(8) = {-4,-12,8,9};

// Side surface, furthest from the z axis
Line Loop(9) = {2,-11,-6,10};

// Top surface
Line Loop(10) = {5,6,7,8};

Plane Surface(20) = {5};
Plane Surface(21) = {6};
Plane Surface(22) = {7};
Plane Surface(23) = {8};
Plane Surface(24) = {9};
Plane Surface(25) = {10};

// Bottom of the box
Physical Line(333) = {5};
// Sides of the box perpendicular to the z axis
Physical Line(666) = {6,7};
// Sides of the box perpendicular to the x axis
Physical Line(777) = {8,9};
// Top of the box
Physical Line(444) = {10};

// This is just to ensure all the interior
// elements get written out. 
Physical Surface(333) = {20};
Physical Surface(666) = {21,22};
Physical Surface(777) = {23,24};
Physical Surface(444) = {25};

Surface Loop(800) = {22, 20, 21, 24, 25, 23};
Volume(900) = {800};
Physical Volume(999) = {900};

