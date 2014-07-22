lx = $DOMAIN_LENGTH_X;
ly = $DOMAIN_LENGTH_Y;
lz = $DOMAIN_LENGTH_Z;
el_sz = $EL_SIZE_X;
sqrt2 = 1.4142135623730951;

Point(1) = {  0,  0,  0, el_sz};
Point(2) = { lx,  0,  0, el_sz};
Point(3) = { lx, ly,  0, el_sz};
Point(4) = {  0, ly,  0, el_sz};
Point(5) = {  0,  0, lz, el_sz};
Point(6) = { lx,  0, lz, el_sz};
Point(7) = { lx, ly, lz, el_sz};
Point(8) = {  0, ly, lz, el_sz};

Point(11) = {  lx/2,  -lx,  0, el_sz};
Point(12) = { lx+ly, ly/2,  0, el_sz};
Point(13) = {  lx/2,    0,  0, el_sz};
Point(14) = {    lx, ly/2,  0, el_sz};
Point(15) = {  lx/2,  -lx, lx, el_sz};
Point(16) = { lx+ly, ly/2, lx, el_sz};
Point(17) = {  lx/2,    0, lx, el_sz};
Point(18) = {    lx, ly/2, lx, el_sz};

Point(21) = {   -lz/sqrt2,    lz/sqrt2, lz/2, el_sz};
Point(22) = { lx-lz/sqrt2,    lz/sqrt2, lz/2, el_sz};
Point(23) = { lx-lz/sqrt2, ly+lz/sqrt2, lz/2, el_sz};
Point(24) = {   -lz/sqrt2, ly+lz/sqrt2, lz/2, el_sz};

Circle(1) = {1, 11, 2};
Circle(2) = {2, 12, 3};
Circle(3) = {3, 13, 4};
Circle(4) = {4, 14, 1};
Circle(5) = {5, 15, 6};
Circle(6) = {6, 16, 7};
Circle(7) = {7, 17, 8};
Circle(8) = {8, 18, 5};
Circle(9) = {1, 21, 5};
Circle(10) = {2, 22, 6};
Circle(11) = {3, 23, 7};
Circle(12) = {4, 24, 8};

Line Loop(21) = {1, 2, 3, 4};
Line Loop(22) = {5, 6, 7, 8};
Line Loop(23) = {10, 6, -11, -2};
Line Loop(24) = {9, -8, -12, 4};
Line Loop(25) = {9, 5, -10, -1};
Line Loop(26) = {12, -7, -11, 3};

Ruled Surface(21) = {21};
Ruled Surface(22) = {22};
Ruled Surface(23) = {23};
Ruled Surface(24) = {24};
Ruled Surface(25) = {25};
Ruled Surface(26) = {26};

Physical Surface(28) = {25};
Physical Surface(29) = {26};
Physical Surface(30) = {21};
Physical Surface(31) = {22};
Physical Surface(32) = {23};
Physical Surface(33) = {24};

Surface Loop(34) = {26, 24, 25, 22, 23, 21};
Volume(1) = {34};
Physical Volume(34) = {1};
