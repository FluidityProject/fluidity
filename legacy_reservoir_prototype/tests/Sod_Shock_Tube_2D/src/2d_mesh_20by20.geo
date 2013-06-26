Point(1) = {0., 0., 0.};
Extrude {1, 0, 0} {
  Point{1}; Layers{20};
}
Extrude {0, 1, 0} {
  Line{1}; Layers{20};
}
Transfinite Surface {5};
// left
Physical Line(6) = {3};
// right
Physical Line(7) = {4};
// bottom
Physical Line(8) = {1};
// top
Physical Line(9) = {2};
// surface
Physical Surface(10) = {5};
