Point(1) = {0,0,0,0.025};
Extrude {1,0,0} {
  Point{1};Layers{50};
}
Extrude {0,1,0} {
  Line{1}; Layers{50};
}
Physical Line(8) = {2};
Physical Line(9) = {4};
Physical Line(10) = {3};
Physical Line(11) = {1};
Physical Surface(7) = {5};

