lx = $DOMAIN_LENGTH_X;
ly = $DOMAIN_LENGTH_Y;
el_sz = $EL_SIZE_X;

Point(1) = {  0,  0, 0, el_sz};
Point(2) = { lx,  0, 0, el_sz};
Point(3) = { lx, ly, 0, el_sz};
Point(4) = {  0, ly, 0, el_sz};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(9) = {1, 2, 3, 4};
Plane Surface(10) = {9};
Physical Line(11) = {1};	// bottom
Physical Line(12) = {2};	// right
Physical Line(13) = {3};	// top
Physical Line(14) = {4};	// left
Physical Surface(15) = {10};

