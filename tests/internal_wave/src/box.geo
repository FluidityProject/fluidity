// Gmsh project created on Thu Mar  6 16:27:33 2008
Point(1) = {0,0,0};
Extrude {0,0,0.6} {
  Point{1}; Layers{30};
}
Extrude {0,0.02,0} {
  Line{1}; Layers{1};
}
Extrude {1.7,0,0} {
  Surface{5}; Layers{80};
}
Physical Surface(38) = {22};
Physical Surface(39) = {14};
Physical Surface(40) = {18};
Physical Surface(41) = {26};
Physical Surface(42) = {5};
Physical Surface(43) = {27};
Physical Volume(44) = {1};
