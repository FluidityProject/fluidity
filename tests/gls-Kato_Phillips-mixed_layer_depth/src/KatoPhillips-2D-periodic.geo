Point(1) = {0,0,0,0.1};
Extrude {250,0,0} {
  Point{1}; Layers{3};
}
Extrude {0,-50,0} {
  Line{1}; Layers{25};
}
Physical Line(8) = {2};
Physical Line(9) = {4};
Physical Line(10) = {3};
Physical Line(11) = {1};
Physical Surface(7) = {5};

