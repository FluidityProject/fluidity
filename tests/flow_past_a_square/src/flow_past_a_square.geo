H = 0.04;
x1 = 5*H + H/2;
x2 = 20*H - H/2;
dx = H/15;

// Outer box
Point(1) = {-x1, -7*H, 0, dx};
Point(2) = {x2, -7*H, 0, dx};
Point(3) = {x2, 7*H, 0, dx};
Point(4) = {-x1, 7*H, 0, dx};

// Inner box
Point(5) = {-H/2, -H/2, 0, dx};
Point(6) = {H/2, -H/2, 0, dx};
Point(7) = {H/2, H/2, 0, dx};
Point(8) = {-H/2, H/2, 0, dx};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {1,2,3,4};

Line(6) = {5,6};
Line(7) = {6,7};
Line(8) = {7,8};
Line(9) = {8,5};
Line Loop(10) = {6,7,8,9};

Plane Surface(11) = {5,10};

// Top at y = h
Physical Line(444) = {3};
// Side along x = -x1 (INLET)
Physical Line(999) = {4};
// Side along x = H+x2 (OUTLET)
Physical Line(666) = {2};
// Bottom at y = 0
Physical Line(333) = {1};

// Top of inner cube at y = H/2
Physical Line(44) = {8};
// Side of inner cube along x = -H/2
Physical Line(99) = {9};
// Side of inner cube along x = H/2
Physical Line(66) = {7};
// Bottom of inner cube at y = -H/2
Physical Line(33) = {6};

// This is just to ensure all the interior
// elements get written out. 
Physical Surface(12) = {11};

