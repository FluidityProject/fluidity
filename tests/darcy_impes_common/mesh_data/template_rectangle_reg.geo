lx = $DOMAIN_LENGTH_X;
ly = $DOMAIN_LENGTH_Y;
nx = $EL_NUM_X;
ny = $EL_NUM_Y;

Point(1) = {0, 0, 0, 0.1};
Extrude {lx, 0, 0} {
  Point{1}; Layers{nx}; 
}
Extrude {0, ly, 0} {
  Line{1}; Layers{ny}; 
}

Physical Line(11) = {1};	// bottom
Physical Line(12) = {4};	// right
Physical Line(13) = {2};	// top
Physical Line(14) = {3};	// left
Physical Surface(15) = {5};

