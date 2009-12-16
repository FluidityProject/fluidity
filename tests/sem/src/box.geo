Point(1) = {0,0,0,1};

Extrude {0,3,0} {
  Point{1};Layers{2};
}

Extrude {0,0,2} {
  Line{1};Layers{6};
}

Line Loop(6) = {2,-4,-1,3};

Plane Surface(7) = {6};

Extrude {0.1,0,0} {
  Surface{5};Layers{1};
}

Physical Surface(55) = {5};
Physical Surface(56) = {16};
Physical Surface(57) = {20,28};
Physical Surface(58) = {24};
Physical Volume (59) = {1};
