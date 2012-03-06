Point (1) = {0, 0, 0, 0.1};
Point (2) = {1, 0, 0, 0.1};
Point (3) = {1, 0, 1.0, 0.1};
Point (4) = {0, 0, 1.0, 0.1};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {3,4,1,2};
Plane Surface(6) = {5};
Extrude {0,0.03,0} {
  Surface{6}; Layers{3};
}
Physical Surface(29) = {15};
Physical Surface(30) = {23};
Physical Surface(31) = {6,28};
Physical Surface(32) = {27,19};
Physical Volume(34) = {1};
