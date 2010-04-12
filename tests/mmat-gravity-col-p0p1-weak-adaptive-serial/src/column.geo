Point (1) = {0, 0, 0};
Extrude {0.5, 0, 0} {
  Point{1}; Layers{10};
}
Extrude {0, 0.5, 0} {
  Line{1}; Layers{10};
}
Extrude {0, 0.5, 0} {
  Line{2}; Layers{10};
}
Physical Line(1) = {1};
Physical Line(2) = {4, 8};
Physical Line(3) = {6};
Physical Line(4) = {7, 3};
Physical Surface(5) = {5};
Physical Surface(6) = {9};
