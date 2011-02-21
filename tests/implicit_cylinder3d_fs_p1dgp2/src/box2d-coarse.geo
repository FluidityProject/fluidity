// box
l1 = 0.1;
l2 = 0.03;

// bottom layer
Point(1) = {-1, -0.31, 0, l1};
Point(2) = {-0.31, -0.31, 0, l2};
Point(3) = {0.31, -0.31, 0, l2};
Point(4) = {1, -0.31, 0, l1};

// top layer
Point(9) = {-1, 0.31, 0, l1};
Point(10) = {-0.31, 0.31, 0, l2};
Point(11) = {0.31, 0.31, 0, l2};
Point(12) = {1, 0.31, 0, l1};

Line(101) = {1, 2};
Line(102) = {2, 3};
Line(103) = {3, 4};
Line(104) = {4, 12};
Line(105) = {12, 11};
Line(106) = {11, 10};
Line(107) = {10, 9};
Line(108) = {9, 1};
Line(109) = {2, 10};
Line(110) = {3, 11};

Line Loop(111) = {107, 108, 101, 109};
Plane Surface(112) = {111};
Line Loop(113) = {106, -109, 102, 110};
Plane Surface(114) = {113};
Line Loop(115) = {105, -110, 103, 104};
Plane Surface(116) = {115};

Physical Line(3) = {108}; // inlet
Physical Line(4) = {104}; // outlet
Physical Line(5) = {107, 106, 105}; // side y+
Physical Line(6) = {101, 102, 103}; // side y-

Physical Surface(117) = {112, 114, 116};
