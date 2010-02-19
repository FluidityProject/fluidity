Point(1) = {0,0,0,1000};
Extrude {1e6,0,0} {
  Point{1}; Layers{10};
}
Extrude {0,1e3,0} {
  Line{1}; Layers{1};
}
Physical Line(1) = {2};
Physical Line(2) = {1};
Physical Line(3) = {3};
Physical Line(4) = {4};
Physical Surface(10) = {5};
