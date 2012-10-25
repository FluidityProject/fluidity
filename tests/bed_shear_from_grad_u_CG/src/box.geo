edgeLength = 0.2;

Point(1) = {0.0, 0.0, 0.0, edgeLength};

Extrude {1, 0.0, 0.0} {
  Point{1}; Layers{1/edgeLength};
}

Extrude {0.0, 1, 0.0} {
  Line{1}; Layers{1/edgeLength};
}

Physical Line(1) = {1,2,3,4};  
Physical Surface(6) = {5};