Point(1) = {0, 0, 0};
Point(2) = {1.0e6, 0, 0};
Point(3) = {1.0e6, 1.0e6, 0};
Point(4) = {0, 1.0e6, 0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {3, 4, 1, 2};
Plane Surface(6) = {5};
Extrude {0, 0, 500} {
  Surface{6};Layers{1};
}
