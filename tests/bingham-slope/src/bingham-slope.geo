Point(1) = {0,0,0,0.1};
Extrude {3.2,0.0,0} {
  Point{1}; Layers{5};
}
Extrude {0.0,1.0,0} {
  Line{1}; Layers{5};
}
Extrude {0.0,0.0,0.1} {
  Surface{5}; Layers{1};
}
Physical Surface(29) = {14};
Physical Surface(30) = {18};
Physical Surface(31) = {22};
Physical Surface(32) = {26};
Physical Surface(33) = {5};
Physical Surface(34) = {27};
Physical Volume(35) = {1};
