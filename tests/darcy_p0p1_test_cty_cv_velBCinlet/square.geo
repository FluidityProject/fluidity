// define a layer variable
lay = 5;

// define a len variable
len = 300.0;

Point(1) = {0.0, 0.0, 0.0, 1.0};

Extrude {len, 0.0, 0.0} {
  Point{1}; Layers{lay}; 
}

Extrude {0.0, len, 0.0} {
  Line{1}; Layers{lay}; 
}

Physical Line(7) = {1};
Physical Line(8) = {4};
Physical Line(9) = {2};
Physical Line(10) = {3};

Physical Surface(11) = {5};
