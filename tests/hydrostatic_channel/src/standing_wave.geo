Point(1) = {0,0,0,1000};
Extrude {1e6,0,0} {
  Point{1}; Layers{10};
}
Extrude {0,1e5,0} {
  Line{1}; Layers{1};
}
Line Loop(6) = {2, -4, -1, 3};
Plane Surface(7) = {6};
//Left end at x=0
Physical Line(3) = {3};
//Right end at x=L
Physical Line(4) = {4};
//Right side at y=0
Physical Line(5) = {1};
//Left side at y=W
Physical Line(6) = {2};
Physical Surface(1) = {5};
