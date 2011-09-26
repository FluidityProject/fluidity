Point(1) = {0, 0, 0};
Extrude {0.01, 0, 0} {
  Point{1}; Layers{50};
}
Extrude {0, 0.01, 0} {
  Line{1}; Layers{50};
}
Physical Surface(1) = {5};
Physical Line(1) = {1};
Physical Line(2) = {3, 4};
Physical Line(3) = {2};