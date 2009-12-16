Point(1) = {0,-10, 0, 0.2};
Point(2) = {0,-10,10, 5.0};

Line(1) = {1,2};

Extrude {20,0,0} {
  Line{1}; Layers{10};
}
Extrude {0,20,0} {
  Surface{5};Layers{10};
}

Physical Surface(28) = {14};
Physical Surface(29) = {26};
Physical Surface(30) = {5,27};
Physical Surface(31) = {18};
Physical Volume(32)  = {1};
