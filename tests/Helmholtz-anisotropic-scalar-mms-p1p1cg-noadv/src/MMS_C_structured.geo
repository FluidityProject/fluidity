Point(1) = {-1.0,-1.0,0};
Extrude {2, 0, 0} {
  Point{1}; Layers{40};
}
Extrude {0, 2, 0} {
  Line{1}; Layers{40};
}
Physical Line(7) = {1};
Physical Line(8) = {4};
Physical Line(9) = {2};
Physical Line(10) = {3};
Physical Surface(11) = {5};

