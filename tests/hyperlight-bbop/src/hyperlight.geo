Point(1) = {0,0,0,1};
Extrude {100,0,0} {
  Point{1};Layers{1};
}
Extrude {0,100,0} {
  Line{1};Layers{1};
}
Extrude {0,0,-100} {
  Surface{5};Layers{50};
}
Physical Surface(111) = {5};
Physical Surface(222) = {14,18,22,26};
Physical Surface(333) = {27};
Physical Volume(444) = {1,2};
Physical Volume(555) = {3};
Physical Volume(666) = {4};
