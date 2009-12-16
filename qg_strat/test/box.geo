Point(1) = {-1,-0.5,0.0,0.1};
Extrude {2,0,0} {
  Point{1}; Layers{31};
}

Extrude {0,1.0,0} {
  Line{1}; Layers{31};
}


Extrude {0,0,0.1} {
  Surface{5}; Layers{10};
}

Physical Surface(1) = {5};
Physical Surface(2) = {27};
Physical Surface(3) = {18,22,26,14};
Physical Volume(1) = {1};
