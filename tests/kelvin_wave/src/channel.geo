Point(1) = {0,0,0,1000};
Extrude {0,100000,0} {
  Point{1}; Layers{10};
}
Extrude {0,0,4000} {
  Line{1}; Layers{1};
}
Extrude {100000,0,0} {
  Surface{5}; Layers{50};
}
Physical Surface(1) = {22};
Physical Surface(2) = {14};
Physical Surface(5) = {5};
Physical Surface(4) = {18,26};
Physical Surface(3) = {27};
Physical Volume(28) = {1};
