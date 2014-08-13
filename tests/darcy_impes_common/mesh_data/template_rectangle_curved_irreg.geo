lx = $DOMAIN_LENGTH_X;
ly = $DOMAIN_LENGTH_Y;
lz = $DOMAIN_LENGTH_Z;
el_sz = $EL_SIZE_X;
sqrt2 = 1.4142135623730951;

Point(1) = {  0,  0,  0, el_sz};
Point(2) = { lx,  0,  0, el_sz};
Point(3) = { lx, ly,  0, el_sz};
Point(4) = {  0, ly,  0, el_sz};
Point(11) = {  lx/2,  -lx,  0, el_sz};
Point(12) = { lx+ly, ly/2,  0, el_sz};
Point(13) = {  lx/2,    0,  0, el_sz};
Point(14) = {    lx, ly/2,  0, el_sz};

Circle(1) = {1, 11, 2};
Circle(2) = {2, 12, 3};
Circle(3) = {3, 13, 4};
Circle(4) = {4, 14, 1};
Line Loop(9) = {1, 2, 3, 4};

Plane Surface(10) = {9};

Physical Line(11) = {1};	// bottom
Physical Line(12) = {2};	// right
Physical Line(13) = {3};	// top
Physical Line(14) = {4};	// left
Physical Surface(15) = {10};
