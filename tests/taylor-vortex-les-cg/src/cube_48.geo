Point(1) = {-1.570796325,-1.570796325,-1.570796325,1};
Extrude {6.2831853, 0, 0} {
  Point{1}; Layers{48};
}
Extrude {0, 6.2831853, 0} {
  Line{1}; Layers{48};
}
Extrude {0, 0, 6.2831853} {
  Surface{5}; Layers{48};
}

//Sides
Physical Surface(1) = {14,22};
//TopBottom
Physical Surface(2) = {5,27};
//Inflow
Physical Surface(3) = {26};
//Outflow
Physical Surface(4) = {18};
//Volume
Physical Volume(1) = {1};
