Point(1) = {0,0,0,0.1};
Extrude {1,0,0} {
  Point{1}; Layers{10};
}
Extrude {0,0,1} {
  Line{1}; Layers{10};
}
Extrude {0,0.05,0.0} {
  Surface{5}; Layers{1};
}
// Front
Physical Surface(28) = {27};
// Back
Physical Surface(29) = {5};
// Left
Physical Surface(30) = {26};
// Top
Physical Surface(31) = {22};
// Right
Physical Surface(32) = {18};
// Bottom
Physical Surface(33) = {14};
// Volume
Physical Volume(34) = {1};
