Point(1) = {0.0, 0.0, 0.0, 1.0};

Extrude {52.6929, 0.0, 0.0} {
  Point{1}; Layers{50}; Recombine;
}

Extrude {0.0, 52.6929, 0.0} {
  Line{1}; Layers{50}; Recombine; 
}

Physical Line(7) = {1};
Physical Line(8) = {4};
Physical Line(9) = {2};
Physical Line(10) = {3};
Physical Surface(10) = {5};
