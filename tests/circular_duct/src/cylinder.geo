l = 0.1;
Point(1) = {1, 0, 0, l};
Point(2) = {0, 0, 0, l};
Point(3) = {-1, 0, 0, l};
Point(4) = {0,1,0,l};
Point(5) = {0,-1,0,l};
Circle(1) = {1, 2, 4};
Circle(2) = {4, 2, 3};
Circle(3) = {3, 2, 5};
Circle(4) = {5, 2, 1};
Line Loop(5) = {3, 4, 1, 2};
Plane Surface(6) = {5};
Extrude {0, 0, 1} {
  Surface{6}; Layers{10};
}
Physical Surface(29) = {6};
Physical Surface(30) = {28};
Physical Surface(31) = {27, 23, 19, 15};
Physical Volume(32) = {1};
