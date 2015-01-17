D = 400; // Diameter of vent
l = 1000; // (Initial) characteristic element length
Lx = 10000;
Ly = 10000;
Lz = 10000;

// Outer frame
Point(1) = {0.0, 0.0, 0.0, l};
Point(2) = {Lx, 0.0, 0.0, l};
Point(3) = {Lx, 0.0, Lz, l};
Point(4) = {0.0, 0.0, Lz, l};
Point(5) = {0.0, Ly, 0.0, l};
Point(6) = {Lx, Ly, 0.0, l};
Point(7) = {Lx, Ly, Lz, l};
Point(8) = {0.0, Ly, Lz, l};

// Inner arc representing the vent
Point(9) = {D/2, 0.0, 0.0, l};
Point(10) = {0.0, 0.0, D/2, l};

Line(1) = {1,9};
Line(2) = {9,2};
Line(3) = {2,3};
Line(4) = {3,4};
Line(5) = {4,10};
Line(6) = {10,1};

Line(7) = {5,6};
Line(8) = {6,7};
Line(9) = {7,8};
Line(10) = {8,5};

Line(11) = {5,1};
Line(12) = {6,2};
Line(13) = {7,3};
Line(14) = {8,4};

Circle(99) = {9, 1, 10};


// Bottom surface
Line Loop(5) = {2,3,4,5,-99};

// Side surface, nearest to the x axis
Line Loop(6) = {1,2,-12,-7,11};

// Side surface, furthest from the x axis
Line Loop(7) = {-4,-13,9,14};

// Side surface, nearest to the z axis
Line Loop(8) = {-6,-5,-14,10,11};

// Side surface, furthest from the z axis
Line Loop(9) = {3,-13,-8,12};

// Top surface
Line Loop(10) = {7,8,9,10};

// Vent
Line Loop(11) = {1,99,6};

Plane Surface(20) = {5};
Plane Surface(21) = {6};
Plane Surface(22) = {7};
Plane Surface(23) = {8};
Plane Surface(24) = {9};
Plane Surface(25) = {10};
Plane Surface(26) = {11};

// Bottom of the box
Physical Line(333) = {2,3,4,5,-99};
// Sides of the box perpendicular to the z axis
Physical Line(666) = {1,2,-12,-7,11,7,-4,-13,9,14};
// Sides of the box perpendicular to the x axis
Physical Line(777) = {-6,-5,-14,10,11,9,3,-13,-8,12};
// Top of the box
Physical Line(444) = {7,8,9,10};
// Vent
Physical Line(999) = {1,99,6};

// This is just to ensure all the interior
// elements get written out. 
Physical Surface(333) = {20};
//Physical Surface(666) = {21,22};
//Physical Surface(777) = {23,24};
Physical Surface(444) = {25};
Physical Surface(999) = {26};

// Side furthest from the x axis
Physical Surface(100) = {22};
// Side nearest to the x axis
Physical Surface(110) = {21};

// Side furthest from the z axis
Physical Surface(200) = {24};
// Side nearest to the z axis
Physical Surface(210) = {23};

Surface Loop(800) = {22, 20, 21, 24, 25, 23, 26};
Volume(900) = {800};
Physical Volume(1000) = {900};

