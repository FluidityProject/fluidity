Point(1) = {0., 0., 0.};
Extrude {4, 0, 0} {
  Point{1}; Layers{10}; 
}
Extrude {0, 0.4, 0} {
  Line{1};  Layers{1};
}
Physical Line(6) = {3};
Physical Line(7) = {4};
Physical Line(8) = {1, 2};
Physical Surface(9) = {5};
