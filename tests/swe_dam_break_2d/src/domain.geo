dx = 5.0;

Point(1) = {0, 0, 0, dx};
Point(2) = {0, 200, 0, dx};
Point(3) = {95, 200, 0, dx};
Point(4) = {95, 170, 0, dx};
Point(5) = {105, 170, 0, dx};
Point(6) = {105, 200, 0, dx};
Point(7) = {200, 200, 0, dx};
Point(8) = {200, 0, 0, dx};
Point(9) = {105, 0, 0, dx};
Point(10) = {105, 95, 0, dx};
Point(11) = {95, 95, 0, dx};
Point(12) = {95, 0, 0, dx};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 1};

Line Loop(5) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};

Plane Surface(888) = {5};

Physical Line(333) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};

Physical Surface(999) = {888};

