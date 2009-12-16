Point(1) = {0,0,0,0.1};
Extrude {1.0,0,0.0} {
  Point{1}; Layers{40}; Recombine; 
}
Extrude {0.0,1.0,0.0} {
  Line{1}; Layers{1}; Recombine; 
}
Physical Line(1) = {3};
Physical Line(2) = {4};
Physical Line(3) = {2,1};
Physical Surface(4) = {5};
