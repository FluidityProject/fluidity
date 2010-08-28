Point(1) = {-1, -1, 0};
Extrude {0, 1, 0} {
  Point{1}; Layers{10};
}
Extrude {0.5, 0, 0} {
  Line{1}; Layers{5};
}
Extrude {0.5, 0, 0} {
  Line{2}; Layers{5};
}
Extrude {0, 1, 0} {
  Line{4}; Layers{10};
}
Extrude {0, 1, 0} {
  Line{8}; Layers{10};
}

// Lower left
Physical Surface(18) = {5};
// Lower right
Physical Surface(19) = {9};
// Upper right
Physical Surface(20) = {17};
// Upper left
Physical Surface(21) = {13};
// Top of the box.
Physical Line(333) = {10, 14};
// Bottom of the box.
Physical Line(444) = {3, 7};
// Side walls.
Physical Line(666) = {6, 16, 11, 1};
