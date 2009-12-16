Point(1) = {0.0,0.0,0,0.1};
Extrude {1,0,0} {
  Point{1};
}
Extrude {0,1,0} {
  Line{1};
}
Extrude {0,0,1} {
  Surface{5};
}
Physical Surface(28) = {27,5,18,26,14,22};
Physical Volume(29) = {1};
