Point(1) = {0, 0, 0};
Extrude {1, 0, 0} {
  Point{1}; Layers{1};
}
Extrude {0, 1, 0} {
  Line{1};  Layers{1};
}
Physical Line(1) = {1};
Physical Line(2) = {4};
Physical Line(3) = {2};
Physical Line(4) = {3};
Physical Surface(5) = {5};
