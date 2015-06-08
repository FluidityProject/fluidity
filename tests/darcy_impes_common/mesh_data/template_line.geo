lx = $DOMAIN_LENGTH_X;
nx = $EL_NUM_X;
dx = $EL_SIZE_X;

Point(1) = {0, 0, 0, dx};

Extrude {lx, 0, 0} {
  Point{1}; Layers{nx}; 
}

// end boundaries
Physical Point(1) = {1};
Physical Point(2) = {2};

Physical Line(3) = {1};

