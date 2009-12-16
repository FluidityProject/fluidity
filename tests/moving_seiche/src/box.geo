// Gmsh project created on Thu Mar  6 16:27:33 2008
Point(1) = {0,0,0};
Extrude {100,0,0} {
  Point{1}; Layers{20};
}
Extrude {0,0,10} {
  Line{1}; Layers{10};
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
Show "*";
Hide {
Volume{1};
}

Show "*";
Hide {
Volume{1};
}

