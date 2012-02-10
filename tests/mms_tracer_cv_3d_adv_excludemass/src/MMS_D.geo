Point(1) = {0.1,-0.3,0.4,0.08};
Extrude {0.5, 0, 0} {
  Point{1}; Layers{48};
}
Extrude {0, 0.4, 0} {
  Line{1}; Layers{40};
}
Extrude {0, 0, 0.3} {
  Surface{5}; Layers{32};
}
Physical Surface(7) = {14};
Physical Surface(8) = {18};
Physical Surface(9) = {22};
Physical Surface(10) = {26};
Physical Surface(11) = {5};
Physical Surface(12) = {27};
Physical Volume(13) = {1};
