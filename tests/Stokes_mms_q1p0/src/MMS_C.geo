pi=3.1415926535897931;
Point(1) = {0., 0., 0., 1.0};
Extrude {pi, 0, 0} {
  Point{1}; Layers{16}; Recombine;
}
Extrude {0, pi, 0} {
  Line{1}; Layers{16}; Recombine;
}
Physical Line(7) = {3};
Physical Line(8) = {4};
Physical Line(9) = {1};
Physical Line(10) = {2};
Physical Surface(11) = {5};
