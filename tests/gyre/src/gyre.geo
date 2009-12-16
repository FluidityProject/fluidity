Point(1) = {0,0,0,0.1};
Extrude {1,0,0} {
  Point{1}; Layers{50};
}
Extrude {0,1,0} {
  Line{1}; Layers{50};
}
Extrude {0,0,0.0005} {
  Surface{5}; Layers{1};
}
Physical Surface(1) = {27};
Physical Surface(2) = {5};
Physical Surface(3) = {14};
Physical Surface(4) = {22};
Physical Surface(5) = {26};
Physical Surface(6) = {18};
Physical Volume(1) = {1};