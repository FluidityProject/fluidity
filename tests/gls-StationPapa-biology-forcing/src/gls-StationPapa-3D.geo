Point(1) = {0,0,0,1};
Extrude {100.0,0.0,0.0} {
  Point{1};Layers{1};
}
Extrude {0.0,100.0,0.0} {
  Line{1};Layers{1};
}
Extrude {0.0,0.0,-150.0} {
  Surface{5};Layers{75};
}
Extrude {0.0,0.0,-48.0} {
  Surface{27};Layers{12};
}

Extrude {0.0,0.0,-44.0} {
  Surface{49};Layers{11};
}
Extrude {0.0,0.0,-50.0} {
  Surface{71};Layers{5};
}
Physical Surface(111) = {5};
// north
Physical Surface(660) = {14,36,58,80};
// south
Physical Surface(670) = {22,44,66,88};
// east
Physical Surface(680) = {18,40,62,84};
// west
Physical Surface(690) = {26,48,70,92};
Physical Surface(333) = {93};
Physical Volume(444) = {1,2};
Physical Volume(555) = {3};
Physical Volume(666) = {4};
