Point(1) = {0, 0, 0, 500};
Point(2) = {13800, 0, 0, 500};
Point(3) = {13800, 1000, 0, 500};
Point(4) = {0, 1000, 0, 500};

Line(1) = {4, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};

Physical Line(20) = {1};
Physical Line(21) = {2};
Physical Line(22) = {3};
Physical Line(23) = {4};

Line Loop(10) = {4, 1, 2, 3};

Plane Surface(11) = {10};
Physical Surface(12) = {11};
