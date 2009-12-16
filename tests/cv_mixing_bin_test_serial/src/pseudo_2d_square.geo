Point(1) = {0.0,0.0,0,0.1};
Point(2) = {3.0,0.0,0,0.1};
Point(3) = {3.0,3.0,0,0.1};
Point(4) = {0.0,3.0,0,0.1};
Line(1) = {4,3};
Line(2) = {3,2};
Line(3) = {2,1};
Line(4) = {1,4};
Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};

Extrude {0,0,0.01} {
  Surface{6}; Layers{1};
}
Physical Surface(29) = {23};
Physical Surface(30) = {19};
Physical Surface(31) = {15};
Physical Surface(32) = {27};
Physical Surface(33) = {28};
Physical Surface(34) = {6};
Physical Volume(35) = {1};
