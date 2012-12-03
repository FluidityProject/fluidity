Point (1) = {0, 0, 0, 0.5};

Extrude {1,0,0} {
  Point{1};Layers{50};Recombine;
}

Extrude {0,1,0} {
  Line{1};Layers{50};Recombine;
}

// Volume number for whole domain.
Physical Surface (1) = {5};
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
