Point (1) = {0, 0, 0, 1.0};
Point (2) = {1.0, 0, 0, 1.0};
Line (1) = {1, 2};

Extrude {0,1.0,0} {
  Line{1}; Layers{1};
}

// Volume number for whole domain.
Physical Surface (1) = {5};
// Free surface
Physical Line(1) = {1};
// Bottom of the box.
Physical Line(2) = {2};
// Side walls.
Physical Line(3) = {3,4};
