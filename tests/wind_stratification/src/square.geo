Point (1) = {0, 0, 0, 1e5};
Point (2) = {1e6, 0, 0, 1e5};
Line (1) = {1, 2};

Extrude {0,-1e6,0} {
  Line{1}; Layers{10};
}

// Volume number for whole domain.
Physical Surface (1) = {5};
// Free surface
Physical Line(1) = {1};
// Bottom of the box.
Physical Line(2) = {2};
// Side walls.
Physical Line(3) = {3,4};
