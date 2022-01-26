radius = 6.37101e+06;
cellSize=1250000;
pio2=Pi/2;

// create inner 1/8 shell
Point(41) = {0, 0, 0, cellSize};
Point(42) = {-radius, 0, 0, cellSize};
Point(43) = {0, radius, 0, cellSize};
// Point(44) = {radius, 0, 0, cellSize};
// Point(45) = {0, -radius, 0, cellSize};
Point(46) = {0, 0, radius, cellSize};
// Point(47) = {0,0,-radius,cellSize};
Circle(21) = {42, 41, 43};
Circle(25)= {46,41,42};
Circle(26)= {46,41,43};
Line Loop(10) = {21, -26, 25} ;
Ruled Surface (60) = {10};

// create outer 1/8 shell
// s0[] = Dilate {{0,0,0}, (radius + cellSize) / radius} {Duplicata{Surface{60};}};

// create remaining 7/8 inner shells
t1[] = Rotate {{0,0,1},{0,0,0},pio2} {Duplicata{Surface{60};}};
t2[] = Rotate {{0,0,1},{0,0,0},pio2*2} {Duplicata{Surface{60};}};
t3[] = Rotate {{0,0,1},{0,0,0},pio2*3} {Duplicata{Surface{60};}};
t4[] = Rotate {{0,1,0},{0,0,0},-pio2} {Duplicata{Surface{60};}};
t5[] = Rotate {{0,0,1},{0,0,0},pio2} {Duplicata{Surface{t4[0]};}};
t6[] = Rotate {{0,0,1},{0,0,0},pio2*2} {Duplicata{Surface{t4[0]};}};
t7[] = Rotate {{0,0,1},{0,0,0},pio2*3} {Duplicata{Surface{t4[0]};}};

// create entire inner and outer shell
Surface Loop(100)={60,t1[0],t2[0],t3[0],t7[0],t4[0],t5[0],t6[0]};
// Surface Loop(200)={s0[0],t1[1],t2[1],t3[1],t7[1],t4[1],t5[1],t6[1]};

// create volume between shells
// Volume(1000)={200,100};

// To create the mesh run
// gmsh sphere.gmsh -2 -v 0 -format msh -o sphere.msh
Physical Surface(101) = {71, 61, 60, 69, 65, 83, 79, 75};
