N=10;
Point (1) = {0, 0, 0};
Extrude {1,0,0} {
  Point{1}; Layers{N};
}
Extrude {0,1,0} {
  Line{1}; Layers{N};
}
Extrude {0,0,1} {
  Surface{5}; Layers{N};
}
Surface Loop(28) = {5, 14, 18, 22, 26, 27};
Volume(29) = {28};
Physical Volume(1) = {29};
