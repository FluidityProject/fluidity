Lx = 0.8;
Ly = 0.1;
Lz = 0.1;

dx = 0.005;

Point(1) = {-Lx, 0.0, 0.0, dx};
Point(2) = {Lx, 0.0, 0.0, dx};
Line(1)={1,2};

Extrude {0, Ly, 0} {
  Line{1};
}

Extrude {0, 0, Lz} {
  Surface{5};
}

// constant y - side walls
Physical Surface(28) = {14, 22};
// constant x - end walls
Physical Surface(29) = {18, 26};
// top
Physical Surface(30) = {27};
//bottom
Physical Surface(31) = {5};

Physical Volume(32) = {1};
