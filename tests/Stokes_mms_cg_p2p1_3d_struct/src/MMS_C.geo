nx = 20;
ny = 20;
nz = 20;
Point(1) = {0, 0, 0};
Extrude {3.1415926535897931, 0, 0} {
  Point{1}; Layers{nx};
}
Extrude {0, 3.1415926535897931, 0} {
  Line{1}; Layers{ny};
}
Extrude {0, 0, 3.1415926535897931} {
  Surface{5}; Layers{nz};
}
Physical Surface(25) = {5};
Physical Surface(26) = {18};
Physical Surface(27) = {27};
Physical Surface(28) = {26};
Physical Surface(29) = {22};
Physical Surface(30) = {14};
Physical Volume(34) = {1};
