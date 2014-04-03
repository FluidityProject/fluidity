dx = 100.0;
Point(1) = {0.0, 0.0, 0.0, dx};

Extrude {dx,0,0} {
  Point{1}; Layers{dx/dx};
}
Extrude {0,10000,0} {
  Line{1}; Layers{10000/dx};
}

// Top
Physical Line(333) = {2};
// Sides
Physical Line(111) = {3};
Physical Line(222) = {4};
// Bottom
Physical Line(999) = {1};

Physical Surface(1000) = {5};
