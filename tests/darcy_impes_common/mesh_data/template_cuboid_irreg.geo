lx = $DOMAIN_LENGTH_X;
ly = $DOMAIN_LENGTH_Y;
lz = $DOMAIN_LENGTH_Z;
el_sz = $EL_SIZE_X;

Point(1) = {  0,  0,  0, el_sz};
Point(2) = { lx,  0,  0, el_sz};
Point(3) = { lx, ly,  0, el_sz};
Point(4) = {  0, ly,  0, el_sz};
Point(5) = {  0,  0, lz, el_sz};
Point(6) = { lx,  0, lz, el_sz};
Point(7) = { lx, ly, lz, el_sz};
Point(8) = {  0, ly, lz, el_sz};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

Line Loop(21) = {1, 2, 3, 4};
Line Loop(22) = {5, 6, 7, 8};
Line Loop(23) = {10, 6, -11, -2};
Line Loop(24) = {9, -8, -12, 4};
Line Loop(25) = {9, 5, -10, -1};
Line Loop(26) = {12, -7, -11, 3};

Plane Surface(21) = {21};
Plane Surface(22) = {22};
Plane Surface(23) = {23};
Plane Surface(24) = {24};
Plane Surface(25) = {25};
Plane Surface(26) = {26};

Physical Surface(28) = {25};
Physical Surface(29) = {26};
Physical Surface(30) = {21};
Physical Surface(31) = {22};
Physical Surface(32) = {23};
Physical Surface(33) = {24};

Surface Loop(34) = {26, 24, 25, 22, 23, 21};
Volume(1) = {34};
Physical Volume(34) = {1};
