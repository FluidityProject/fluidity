
Point(1) = {0, 0, 0, 0.0625};
Extrude {-1, 1, 0} {
  Point{1};Layers{16};
}
Point(3) = {1, 1, 0, 0.0625};
Extrude {-1, 1, 0} {
  Point{3};Layers{16};
}
Line(3)={1,3};
Line(4)={2,4};
Line Loop(5) = {4, -2, -3, 1};
Plane Surface(6) = {5};

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {4, 3};
Physical Surface(1) = {6};
