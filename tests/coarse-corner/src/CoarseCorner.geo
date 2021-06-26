Point(1) = {0.0,0.0,0.0,0.1};
Extrude {10.0,0,0} {
  Point{1}; Layers{10};
}Extrude {0,10.0,0} {
  Line{1}; Layers{10};
}
Extrude {0,0,10} {
  Surface{5}; Layers{10};
}
Physical Volume(28) = {1};
Physical Surface(29) = {18,5,22,26,27,14};
