Point (1) = {-3e6, -3e6, 0, 1e5};

Extrude {6e6,0,0} {
  Point{1};Layers{60};Recombine;
}

Extrude {0,6e6,0} {
  Line{1};Layers{60};Recombine;
}

// Volume number for whole domain.
Physical Surface (1) = {5};
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
