Point(1) = {0,0,0,0.05};
Extrude {20,0,0} {
  Point{1};
}
Extrude {0,1,0} {
  Line{1};
}
Physical Line(333) = {1,3,2,4};
Physical Surface(1) = {5};
