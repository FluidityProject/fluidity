s=0.02;
b=0.1;

Point(1) = {0.00, 0.0, 0.0, b};
Point(2) = {0.25, 0.0, 0.0, s};
Point(3) = {1.00, 0.0, 0.0, b};

Point(4) = {0.00, 0.5, 0.0, s};
Point(5) = {0.25, 0.5, 0.0, s};
Point(6) = {1.00, 0.5, 0.0, b};

Point(7) = {0.00, 1.0, 0.0, b};
Point(8) = {0.25, 1.0, 0.0, b};
Point(9) = {1.00, 1.0, 0.0, b};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {4, 5};
Line(4) = {5, 6};
Line(5) = {7, 8};
Line(6) = {8, 9};
Line(7) = {7, 4};
Line(8) = {4, 1};
Line(9) = {9, 6};
Line(10) = {6, 3};
Line(11) = {2, 5};
Line(12) = {5, 8};
Line Loop(13) = {8, 1, 11, -3};
Plane Surface(14) = {13};
Line Loop(15) = {2, -10, -4, -11};
Plane Surface(16) = {15};
Line Loop(17) = {4, -9, -6, -12};
Plane Surface(18) = {17};
Line Loop(19) = {5, -12, -3, -7};
Plane Surface(20) = {19};

// left
Physical Line(6) = {7, 8};
// right
Physical Line(7) = {10, 9};
// bottom
Physical Line(8) = {1, 2};
// top
Physical Line(9) = {5, 6};
// surface
Physical Surface(10) = {14, 16, 18, 20};
