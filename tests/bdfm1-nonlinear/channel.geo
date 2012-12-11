Ls=50;
Point(1) = {0, 0, 0, 0.1};
Extrude {6.283185307179586, 0, 0} {
  Point{1};Layers{Ls};
}
Extrude {0., 6.283185307179586, 0} {
  Line{1}; Layers{Ls};
}
Physical Line(1) = {3};
Physical Line(2) = {4};
Physical Line(3) = {1};
Physical Line(4) = {2};
Physical Surface(5) = {5};
