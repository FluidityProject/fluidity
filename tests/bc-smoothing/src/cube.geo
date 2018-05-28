Point(1) = {0,0,0,1.0};
Extrude {1,0,0} {
  Point{1}; Layers{10};
}
Extrude {0,1,0} {
  Line{1}; Layers{10};
}
Extrude {0,0,1} {
  Surface{5}; Layers{10};
}
Physical Volume(1) = {1};
// Left x=0
Physical Surface(1) = {26};
// Right x=1
Physical Surface(2) = {18};
// Front y=0
Physical Surface(3) = {14};
// Back y=1
Physical Surface(4) = {22};
// Bottom z=0
Physical Surface(5) = {5};
// Top z=1
Physical Surface(6) = {27};
