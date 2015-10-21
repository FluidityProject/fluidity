lto=1; // number of cells in a pi/4 segment in the tangential direction of the inner radius
p=1;
Point(1) = {1.2, 0, 0, 1.0};

Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{1};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{2};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{4};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{5};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{6};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{7};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{8};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{9};
}
Transfinite Line {1} = lto Using Progression p;
Transfinite Line {2} = lto Using Progression p;
Transfinite Line {3} = lto Using Progression p;
Transfinite Line {4} = lto Using Progression p;
Transfinite Line {5} = lto Using Progression p;
Transfinite Line {6} = lto Using Progression p;
Transfinite Line {7} = lto Using Progression p;
Transfinite Line {8} = lto Using Progression p;
Physical Line(9) = {1};
Physical Line(10) = {2};
Physical Line(11) = {3};
Physical Line(12) = {4};
Physical Line(13) = {5};
Physical Line(14) = {6};
Physical Line(15) = {7};
Physical Line(16) = {8};
