lx = $DOMAIN_LENGTH_X;
ly = $DOMAIN_LENGTH_Y;
lz = $DOMAIN_LENGTH_Z;
nx = $EL_NUM_X;
ny = $EL_NUM_Y;
nz = $EL_NUM_Z;

Point(1) = {0, 0, 0, 0.1};
Extrude {lx, 0, 0} {
  Point{1}; Layers{nx}; 
}
Extrude {0, 0, lz} {
  Line{1}; Layers{nz}; 
}
Extrude {0, ly, 0} {
  Surface{5}; Layers{ny};
}

Physical Surface(28) = {5};
Physical Surface(29) = {27};
Physical Surface(30) = {26};
Physical Surface(31) = {22};
Physical Surface(32) = {18};
Physical Surface(33) = {14};
Physical Volume(34) = {1};

