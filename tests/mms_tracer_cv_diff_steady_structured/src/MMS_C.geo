Point(1) = {0.1,-0.3,0,0.1};
Extrude {0.5,0,0.0} {
  Point{1}; Layers{41}; 
}
Extrude {0.0,0.4,0.0} {
  Line{1}; Layers{33}; 
}
Physical Line(7) = {1};
Physical Line(8) = {4};
Physical Line(9) = {2};
Physical Line(10) = {3};
Physical Surface(11) = {5};
