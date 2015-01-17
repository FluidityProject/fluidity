D = 500.0;
Lx = 40000.0;
Ly = 10000.0;
dx = 100;

Point(1) = {0.0, 0.0, 0.0, dx};
Point(2) = {0.0, Ly, 0.0, dx};
Point(3) = {Lx, Ly, 0.0, dx};
Point(4) = {Lx, 0.0, 0.0, dx};
Point(5) = {D/2, 0.0, 0.0, dx};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,1};
Line Loop(6) = {1,2,3,4,5};

Plane Surface(7) = {6};

// Left side
Physical Line(111) = {1};
// Right side
Physical Line(222) = {3};
// Top
Physical Line(333) = {2};
// Bottom (EXCLUDING VENT)
Physical Line(444) = {4};
// Bottom (VENT)
Physical Line(999) = {5};

// This is just to ensure all the interior
// elements get written out. 
Physical Surface(8) = {7};
