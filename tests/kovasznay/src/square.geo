dx = 0.05;
Point (1) = {-0.5, -0.5, 0, dx};
Point (2) = {1, -0.5, 0, dx};
Point (3) = {1, 1.5, 0, dx};
Point (4) = {-0.5, 1.5, 0, dx};
Line(1) = {4, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

Physical Line(7) = {1};
Physical Line(8) = {2};
Physical Line(9) = {3};
Physical Line(10) = {4};
Physical Surface(11) = {6};
