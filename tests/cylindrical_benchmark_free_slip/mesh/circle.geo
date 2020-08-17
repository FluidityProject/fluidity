lti = 16*N; // number of cells per Pi/4 section
Point(1) = {2.22, 0, 0, 1.0};
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
Transfinite Line {1} = lti + 1;
Transfinite Line {2} = lti + 1;
Transfinite Line {3} = lti + 1;
Transfinite Line {4} = lti + 1;
Transfinite Line {5} = lti + 1;
Transfinite Line {6} = lti + 1;
Transfinite Line {7} = lti + 1;
Transfinite Line {8} = lti + 1;
Physical Line(9) = {8, 7, 6, 5, 4, 3, 2, 1};
