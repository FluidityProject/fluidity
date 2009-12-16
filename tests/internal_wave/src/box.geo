// Gmsh project created on Thu Mar  6 16:27:33 2008
Point(1) = {0,0,0};
Extrude {1.7,0,0} {
  Point{1}; Layers{80};
}
Extrude {0,0,0.6} {
  Line{1}; Layers{30};
}
Extrude {0,0.02,0} {
  Surface{5}; Layers{1};
}
Physical Surface(38) = {5};
Physical Surface(39) = {27};
Physical Surface(40) = {22};
Physical Surface(41) = {14};
Physical Surface(42) = {26};
Physical Surface(43) = {18};
Physical Volume(44) = {1};
