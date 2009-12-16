Point(1) = {0,0,0,0.1};
Extrude {0,0,1} {
  Point{1}; Layers{5};
}
Extrude {0,2,0} {
  Line{1}; Layers{5};
}
Extrude {4,0,0} {
  Surface{5}; Layers{5};
}
Physical Surface(1) = {5};
Physical Surface(2) = {27};
Physical Surface(3) = {18,14,26,22};

Physical Volume(31) = {1};
