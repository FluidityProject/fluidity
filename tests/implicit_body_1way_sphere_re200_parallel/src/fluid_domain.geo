// Element lengths for the zones of meshes

//sim8_3 (Re=100):

// lenghts of elements of boxes/zones in x-direction:
lc = 7;
lm = 0.5;
lf = 0.09;
lc2 = 7;

// variables for the size/height of boxes/zones in y- and z-direction:
den = 2*2;
yc = 4 /den;
ym = 2 /den;
yf = 1.5 /den;
zc = yc;
zm = ym;
zf = yf;
// variables for the layers in y- and z-direction of the boxes/zones:
//laf = 20;
//lac = yc*10/lc;
lac = 3;
//lam = ym*10/lm;
lam = 4;
//laf = yf*10/lf;
laf = 20;


// Center line (with different element lengths)
Point(1) = {-2,0,0,lc};
Point(2) = {-1,0,0,lm};
Point(3) = {-0.75,0,0,lf};
Point(4) = { 0.75,0,0,lf};
Point(5) = { 1,0,0,lm};
Point(6) = {2,0,0,lc2};
// Center line:
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {5,4};
Line(5) = {6,5};



// Extruding in Y:
Extrude {0,2,0} {
  Line{1,2,3,4,5}; Layers{ {laf,lam,lac} , {yf, ym, yc}}; 
}
Extrude {0,-2,0} {
  Line{1,2,3,4,5}; Layers{ {laf,lam,lac} , {yf, ym, yc}};
}

// Extruding in Z:
Extrude {0,0,2} {
  Surface{9, 13, 17, 21, 25}; Layers{ {laf,lam,lac} , {zf, zm, zc}};
}
Extrude {0,0,-2} {
  Surface{9, 13, 17, 21, 25}; Layers{ {laf,lam,lac} , {zf, zm, zc}};
}
Extrude {0,0,2} {
  Surface{29, 33, 37, 41, 45}; Layers{ {laf,lam,lac} , {zf, zm, zc}};
}
Extrude {0,0,-2} {
  Surface{29, 33, 37, 41, 45}; Layers{ {laf,lam,lac} , {zf, zm, zc}};
}


// Physical IDs:
//front/inlet:
Physical Surface(1) = {66, 176, 286, 396};
//sides:
Physical Surface(2) = {172, 62, 84, 194, 216, 106, 238, 128, 260, 150, 282, 392, 304, 414, 436, 326, 458, 348, 370, 480};
//top and bottom
Physical Surface(3) = {67, 287, 89, 309, 331, 111, 133, 353, 375, 155, 177, 397, 419, 199, 221, 441, 463, 243, 265, 485};
//back/outlet:
Physical Surface(4) = {154, 264, 484, 374};

//volumes
Physical Volume(5) = {1, 11, 16, 6, 7, 17, 12, 2, 3, 13, 18, 8, 9, 19, 14, 4, 5, 15, 20, 10};

