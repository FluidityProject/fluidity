Point (1) = {0.0, 0.25, 0.0, 0.025};
Extrude {0.125,0,0} {
  Point{1}; Layers{1}; Recombine;
}
Extrude {0,0.5,0} {
  Line{1}; Layers{10}; Recombine;
}
Physical Line(6) = {1};
Physical Line(7) = {4};
Physical Line(8) = {2};
Physical Line(9) = {3};
Physical Surface(10) = {5};
