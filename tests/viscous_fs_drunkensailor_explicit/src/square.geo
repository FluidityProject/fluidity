Point(1) = {-250e3, -500e3, 0, 2e4};
Extrude {500.e3, 0, 0} {
  Point{1}; Layers{50};
}
Extrude {0.0, 500.e3, 0} {
  Line{1}; Layers{50};
}
Physical Line(6) = {1};
Physical Line(7) = {4};
Physical Line(8) = {2};
Physical Line(9) = {3};
Physical Surface(10) = {5};
