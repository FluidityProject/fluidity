Point(1) = {0., 0., 0., 0.1};
Extrude {0.1, 0, 0} {
  Point{1}; Layers{1};
}
Extrude {0, 0.1, 0} {
  Line{1}; Layers{1};
}

Extrude {5, 0, 0} {
  Line{4}; Layers{10};
}
Extrude {0, 5, 0} {
  Line{2}; Layers{10};
}
Extrude {5, 0, 0} {
  Line{12}; Layers{10};
}
// short left
Physical Line(18) = {3};
// short bottom
Physical Line(19) = {1};
// right
Physical Line(20) = {6, 14};
// top
Physical Line(21) = {16, 10};
// long left
Physical Line(22) = {11};
// long bottom
Physical Line(23) = {7};
// surface
Physical Surface(24) = {17, 13, 9, 5};
