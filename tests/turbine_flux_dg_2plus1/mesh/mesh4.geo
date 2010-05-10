Point(1) = {0, 0, 0, 0.1};
Point(2) = {1, 0, 0, 0.1};
Point(3) = {1, 1, 0, 0.1};
Point(4) = {0, 1, 0, 0.1};
Point(5) = {0.5, 0, 0, 0.1};
Point(6) = {0.5, 1, 0, 0.1};
Point(7) = {0.5000000001, 0, 0, 0.1};
Point(8) = {0.5000000001, 1, 0, 0.1};

Line(1) = {4, 1};
Line(2) = {1, 5};
Line(3) = {5, 6};
Line(4) = {6, 4};

Line(5) = {8, 7};
Line(6) = {7, 2};
Line(7) = {2, 3};
Line(8) = {3, 8};

Physical Line(20) = {1};
Physical Line(21) = {2};
Physical Line(24) = {3};
Physical Line(22) = {4};

Physical Line(25) = {5};
Physical Line(26) = {6};
Physical Line(27) = {7};
Physical Line(28) = {8};

Line Loop(10) = {1, 2, 3, 4};
Line Loop(11) = {5, 6, 7, 8};

Plane Surface(12) = {10};
Plane Surface(13) = {11};
Physical Surface(14) = {12, 13};
