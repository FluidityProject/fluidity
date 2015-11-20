Point(1) = {0, 0, 0, 0.1};
Extrude {1, 0, 0} {
  Point{1}; Layers{10};
}
Extrude {0, 0.5, 0} {
  Line{1}; Layers{5};
}
Extrude {0, 0.5, 0} {
  Line{2}; Layers{5};
}
Physical Line(1) = {1, 4, 8, 6, 7, 3};
Physical Surface(1) = {5}; // bottom half
Physical Surface(2) = {9}; // top half
