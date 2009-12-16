Point(1) = {0, 0, 0};
Extrude {0, 1, 0} {
  Point{1};Layers{1};
}
Extrude {0, 0, 1} {
  Line{1};Layers{1};
}
Extrude {1, 0, 0} {
  Surface{5};Layers{5};
}
Physical Surface(1) = {5};
Physical Surface(2) = {27};
Physical Surface(3) = {26};
Physical Surface(4) = {18};
Physical Surface(5) = {22};
Physical Surface(6) = {14};
Physical Volume(28) = {1};
