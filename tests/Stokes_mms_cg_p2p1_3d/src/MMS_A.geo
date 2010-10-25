nx = 5;
ny = 5;
nz = 5;
Point(1) = {0, 0, 0};
Extrude {3.1415926535897931, 0, 0} {
  Point{1}; Layers{nx};
}
Extrude {0, 3.1415926535897931, 0} {
  Point{2}; Layers{ny};
}
Extrude {-3.1415926535897931, 0, 0} {
  Point{3}; Layers{nx};
}
Extrude {0, -3.1415926535897931, 0} {
  Point{4}; Layers{ny};
}
Extrude {0, 0, 3.1415926535897931} {
  Point{1}; Layers{nz};
}
Extrude {0, 0, 3.1415926535897931} {
  Point{2}; Layers{nz};
}
Extrude {0, 0, 3.1415926535897931} {
  Point{3}; Layers{nz};
}
Extrude {0, 0, 3.1415926535897931} {
  Point{4}; Layers{nz};
}
Extrude {3.1415926535897931, 0, 0} {
  Point{5}; Layers{nx};
}
Extrude {0, 3.1415926535897931, 0} {
  Point{6}; Layers{ny};
}
Extrude {-3.1415926535897931, 0, 0} {
  Point{7}; Layers{nx};
}
Extrude {0, -3.1415926535897931, 0} {
  Point{8}; Layers{ny};
}
Line Loop(13) = {1, 2, 3, 4};
Plane Surface(14) = {13};
Line Loop(15) = {6, 10, -7, -2};
Plane Surface(16) = {15};
Line Loop(17) = {10, 11, 12, 9};
Plane Surface(18) = {17};
Line Loop(19) = {12, -5, -4, 8};
Plane Surface(20) = {19};
Line Loop(21) = {3, 8, -11, -7};
Plane Surface(22) = {21};
Line Loop(23) = {1, 6, -9, -5};
Plane Surface(24) = {23};
Physical Surface(25) = {14};
Physical Surface(26) = {16};
Physical Surface(27) = {18};
Physical Surface(28) = {20};
Physical Surface(29) = {22};
Physical Surface(30) = {24};
Surface Loop(31) = {24, 14, 16, 18, 22, 20};
Volume(32) = {31};
Physical Volume(33) = {32};
