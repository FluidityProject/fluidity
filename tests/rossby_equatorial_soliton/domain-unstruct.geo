Point(1) = {-24.0,-8.0,0.0,0.15};
Point(2) = {-24.0,8.0,0.0,0.15};
Point(3) = {24.0,8.0,0.0,0.15};
Point(4) = {24.0,-8.0,0.0,0.15};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(5) = {3,4,1,2};
Plane Surface(6) = {5};
Extrude {0,0,1} {
  Surface{6}; Layers{1};
}


Physical Surface(28) = {6};
Physical Surface(29) = {28};
Physical Surface(30) = {27,19};
Physical Surface(31) = {23};
Physical Surface(32) = {15};
Physical Volume(34) = {1};
