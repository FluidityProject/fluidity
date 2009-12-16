Point(1) = {0,0,0,0.5};
Extrude {1,1,0} {
  Point{1};
}
Extrude {-1,1,0} {
  Line{1};
}
Physical Line(6) = {4,2,3,1};
Physical Surface(7) = {5};
