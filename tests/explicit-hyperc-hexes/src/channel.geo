Point (1) = {0, -0.5, 0, 0.025};
Extrude {0.2,0,0} {
  Point{1}; Layers{1}; Recombine;
}
Extrude {0,1.0,0} {
  Line{1}; Layers{20}; Recombine;
}
Extrude {0,0,0.2} {
  Surface{5}; Layers{1}; Recombine;
}
Physical Surface(29) = {14};
Physical Surface(30) = {18};
Physical Surface(31) = {22};
Physical Surface(32) = {26};
Physical Surface(33) = {27};
Physical Surface(34) = {5};
Physical Volume(35) = {1};
