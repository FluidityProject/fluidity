Point(1) = {-0.5,-0.5,0,0.05};
Extrude {1,0,0} {
  Point{1};
}
Extrude {0,2.5,0} {
  Line{1};
}
Extrude {0,0,0.1} {
  Surface{5};
}
Physical Surface(6) = {14};
Physical Surface(7) = {18};
Physical Surface(8) = {22};
Physical Surface(9) = {26};
Physical Surface(10) = {5};
Physical Surface(11) = {27};
Physical Volume(12) = {1};
