N=20;
Point (1) = {0, 0, 0};
Extrude {1,0,0} {
  Point{1}; Layers{N};
}
Extrude {0,1,0} {
  Line{1}; Layers{N};
}
Line Loop(6) = {2, -4, -1, 3};
Plane Surface(7) = {6};
Physical Surface(1) = {7};
