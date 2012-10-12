Point(1) = {0., 0., 0., 0.03};
Extrude {1, 0, 0} {
  Point{1};
}
Extrude {0, 0.1, 0} {
  Line{1};
}
// left
Physical Line(6) = {3};
// right
Physical Line(7) = {4};
// top and bottom
Physical Line(8) = {1, 2};
// volume
Physical Surface(9) = {5};
