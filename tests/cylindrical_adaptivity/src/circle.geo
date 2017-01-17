n=2; // number of cells in tangential direction in each Pi/4 section
Point(1) = {2.22, 0, 0, 1.0};
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{1}; Layers{n};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{2}; Layers{n};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{4}; Layers{n};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{5}; Layers{n};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{6}; Layers{n};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{7}; Layers{n};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{8}; Layers{n};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Point{9}; Layers{n};
}
Delete {
  Point{3};
}
Physical Line(9) = {1, 2, 3, 4, 5, 6, 7, 8};
