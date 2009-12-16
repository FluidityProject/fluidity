Point(1) = {0,0,0,1000};
Extrude {1e6,0,0} {
  Point{1}; Layers{10};
}
Extrude {0,1e5,0} {
  Line{1}; Layers{1};
}
Extrude {0,0,1e3} {
  Surface{5}; Layers{1};
}
Physical Surface(1) = {27};
Physical Surface(2) = {5};
Physical Surface(3) = {18};
Physical Surface(4) = {26};
Physical Surface(5) = {22};
Physical Surface(6) = {14};
Physical Volume(1) = {1};