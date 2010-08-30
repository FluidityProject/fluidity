layers=5;
Point(1) = {0, 0, 0, 0.25};
Extrude {0.5, 0, 0} {
  Point{1}; Layers{layers};
}
Extrude {0.5, 0, 0} {
  Point{2}; Layers{layers};
}
Extrude {0, 0.5, 0} {
  Point{3}; Layers{layers};
}
Extrude {0, 0.5, 0} {
  Point{4}; Layers{layers};
}
Extrude {0, 0.5, 0} {
  Point{1}; Layers{layers};
}
Extrude {0, 0.5, 0} {
  Point{6}; Layers{layers};
}
Extrude {0.5, 0, 0} {
  Point{7}; Layers{layers};
}
Extrude {0.5, 0, 0} {
  Point{8}; Layers{layers};
}
Extrude {0.5, 0, 0} {
  Point{6}; Layers{layers};
}
Extrude {0.5, 0, 0} {
  Point{9}; Layers{layers};
}
Extrude {0, 0.5, 0} {
  Point{2}; Layers{layers};
}
Extrude {0, 0.5, 0} {
  Point{9}; Layers{layers};
}
Line Loop(13) = {10, -3, -2, 11};
Plane Surface(14) = {13};
Line Loop(15) = {11, -9, -5, 1};
Plane Surface(16) = {15};
Line Loop(17) = {9, 12, -7, -6};
Plane Surface(18) = {17};
Line Loop(19) = {12, 8, -4, -10};
Plane Surface(20) = {19};
Physical Line(21) = {1, 2};
Physical Line(22) = {3, 4};
Physical Line(23) = {8, 7};
Physical Line(24) = {6, 5};
Physical Surface(25) = {16};
Physical Surface(26) = {14};
Physical Surface(27) = {20};
Physical Surface(28) = {18};
