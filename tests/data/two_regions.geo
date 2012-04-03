Point(1) = {0,0,0,0.1};
Extrude {1,0,0} {
  Point{1}; Layers{10};
}
Extrude {0,1,0} {
  Line{1}; Layers{10};
}
Extrude {0,1,0} {
  Line{2}; Layers{10};
}
Physical Line(10) = {6, 8, 4, 1, 3, 7};
Physical Surface(11) = {9};
Physical Surface(12) = {5};
