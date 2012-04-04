Point(1) = {0., 0., 0., 0.03};
Extrude {1, 0, 0} {
  Point{1}; 
}
Extrude {0, 0.1, 0} {
  Line{1}; 
}
Physical Line(6) = {3};
Physical Line(7) = {4};
Physical Line(8) = {1, 2};
Physical Surface(9) = {5};
