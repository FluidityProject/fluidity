Point(1) = {0., 0., 0., 0.1};
Extrude {1, 0, 0} {
  Point{1}; Layers{4};
}
Extrude {0, 0.1, 0} {
  Line{1}; Layers{2};
}
Extrude {0, 0, 0.1} {
  Surface{5}; Layers{2};
}
// Inflow
Physical Surface(1) = {26};
// Sides
Physical Surface(2) = {14, 22};
// Top-Bottom
Physical Surface(3) = {5, 27};
// Outflow
Physical Surface(4) = {18};
// Volume
Physical Volume(5) = {1};
