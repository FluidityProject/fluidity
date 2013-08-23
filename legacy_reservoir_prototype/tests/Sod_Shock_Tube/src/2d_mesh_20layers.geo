Point(1) = {0., 0., 0.};
Extrude {1, 0, 0} {
  Point{1}; Layers{40};
}
Extrude {0, 0.025, 0} {
  Line{1}; Layers{1};
}
Physical Line(6) = {3};
Physical Line(7) = {4};
Physical Line(8) = {1, 2};
Physical Surface(9) = {5};
