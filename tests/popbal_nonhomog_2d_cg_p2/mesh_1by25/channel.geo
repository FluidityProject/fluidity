edgeLength = 1.0/25.0;

Point(1) = {0.0, 0.0, 0.0, 1};

Extrude {1.0, 0.0, 0.0} {
  Point{1}; Layers{1.0/edgeLength};
}

Extrude {0.0, 1.0, 0.0} {
  Line{1}; Layers{1.0/edgeLength};
}

Physical Line(1) = {1}; 
Physical Line(2) = {2};  
Physical Line(3) = {3};  
Physical Line(4) = {4};  
Physical Surface(6) = {5};
