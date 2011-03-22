// define a layer variable
lay = 10;

// define a len variable
len = 100.0;

Point (1) = {0.0, 0.0, 0.0, 1.0};

Extrude {len, 0.0, 0.0} {
  Point{1}; Layers{lay}; Recombine;
}

Extrude {0.0, len, 0.0} {
  Line{1}; Layers{lay}; Recombine;
}

Extrude {0.0, 0.0, len} {
  Surface{5}; Layers{lay}; Recombine;
}

Physical Surface(29) = {14};
Physical Surface(30) = {18};
Physical Surface(31) = {22};
Physical Surface(32) = {26};
Physical Surface(33) = {27};
Physical Surface(34) = {5};

Physical Volume(35) = {1};
