Point(1) = {0, 0, 0};
Extrude {1, 0, 0} {
  Point{1};
}
Extrude {0, 1, 0} {
  Line{1};
}
// Left
Physical Line(1) = {3};
// Right
Physical Line(2) = {4};
// Bottom
Physical Line(3) = {1};
// Top
Physical Line(4) = {2};
Physical Surface(1) = {5};
