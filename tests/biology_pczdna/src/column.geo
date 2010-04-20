Point(1) = {0,0,0,1};
Extrude {100,0,0} {
  Point{1};Layers{1};
}
Extrude {0,100,0} {
  Line{1};Layers{1};
}
Extrude {0,0,-200} {
  Surface{5};Layers{200};
}
Extrude {0,0,-150} {
  Surface{27};Layers{150};
}

Extrude {0,0,-50} {
  Surface{49};Layers{50};
}
Extrude {0,0,-100} {
  Surface{71};Layers{100};
}
Physical Surface(111) = {5};
Physical Surface(222) = {18,22,14,26,36,48,40,44,58,70,62,66,80,84,88,92};
Physical Surface(333) = {93};
Physical Volume(444) = {1,2};
Physical Volume(555) = {3};
Physical Volume(666) = {4};
