pi = 3.141592653589793;
n = XX;

Point(1) = {0.0, 0.0, 0.0, 1};

Extrude {pi, 0.0, 0.0} {
  Point{1}; Layers{n};
}

Extrude {0.0, pi, 0.0} {
  Line{1}; Layers{n};
}

Physical Line(1) = {1,2,3,4}; 
Physical Surface(6) = {5};
