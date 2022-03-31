Point(1) = {-0.1,0.2,0,0.1};
Extrude {0.8,0,0.0} {
  Point{1}; Layers{11}; Recombine;
}
Extrude {0.0,0.6,0.0} {
  Line{1}; Layers{9}; Recombine;
}
Physical Line(7) = {1};
Physical Line(8) = {4};
Physical Line(9) = {2};
Physical Line(10) = {3};
Physical Surface(11) = {5};
