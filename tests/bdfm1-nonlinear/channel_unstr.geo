
Point(1) = {0, 0, 0, 6.283185307179586/50};
Extrude {0, 6.283185307179586, 0} {
  Point{1};Layers{100};
}
Point(3) = {6.283185307179586, 0, 0, 6.283185307179586/50};
Extrude {0, 6.283185307179586, 0} {
  Point{3};Layers{100};
}
Line(3)={1,3};
Line(4)={2,4};
Line Loop(5) = {4, -2, -3, 1};
Plane Surface(6) = {5};

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
Physical Surface(1) = {6};

