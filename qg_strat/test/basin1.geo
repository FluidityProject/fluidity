Point(1) = {-500000,-500000,0.0,200000.0};
Extrude {1000000.0,0,0} {
  Point{1}; Layers{31};
}

Extrude {0,1000000.0,0} {
  Line{1}; Layers{31};
}


Extrude {0,0,1000.0} {
  Surface{5}; Layers{10};
}

Physical Surface(1) = {5};
Physical Surface(2) = {27};
Physical Surface(3) = {18,22,26,14};
Physical Volume(1) = {1};
