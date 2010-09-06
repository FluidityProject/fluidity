Point(1) = {0, 0, 0};
Extrude {1, 0, 0} {
  Point{1};
}
Extrude {0, 1, 0} {
  Line{1};
}
Extrude {0, 0, 1} {
  Surface{5};
}
Physical Surface(28) = {26};  // Left
Physical Surface(29) = {18};  // Right
Physical Surface(30) = {14};  // Front
Physical Surface(31) = {22};  // Back
Physical Surface(32) = {5};   // Bottom
Physical Surface(33) = {27};  // Top
Physical Volume(34) = {1};
