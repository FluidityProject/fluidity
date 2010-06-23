dx = 0.1;
Point (1) = {0.0, 0.0, 0, dx};
Point (2) = {1.0, 0.0, 0, dx};
Point (3) = {0.5, Sqrt(3.0/4.0), 0, dx};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 1};
Line Loop(4) = {3, 1, 2};
Plane Surface(5) = {4};
Extrude {0, 0, 1} {
  Surface{5};
}
Line Loop(23) = {12, -7, -11, 3};
Surface Loop(24) = {17, 5, 13, 22, 21};
Volume(25) = {24};
Physical Surface(26) = {5};
Physical Surface(27) = {22};
Physical Surface(28) = {17};
Physical Surface(29) = {21};
Physical Surface(30) = {13};
Physical Volume(31) = {1};
