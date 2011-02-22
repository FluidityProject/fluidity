lc= 0.03;

Point(1) = {0,0,0.0,lc};
Point(2) = {0.14,0,0.0,lc};
Point(3) = {0,0.14,0.0,lc};
Point(4) = {-0.14,0,0.0,lc};
Point(5) = {0,-0.14,0.0,lc};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

Extrude {0, 0, -0.5} {
  Surface{6};
Layers{5};
}

Physical Surface(29) = {6, 19, 23, 27, 15, 28};
Physical Volume(30) = {1};
