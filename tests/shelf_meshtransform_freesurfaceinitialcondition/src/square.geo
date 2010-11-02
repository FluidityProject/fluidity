Point (1) = {0, 0, 0, 5e3};
Point (2) = {5e5, 0, 0, 5e4};
Point (3) = {1e6, 0, 0, 5e4};
Line (1) = {1, 2};
Line (2) = {2, 3};

Extrude {0,-1e6,0} {
  Line{1,2}; Layers{1};
}

// Volume number for whole domain.
Physical Surface (1) = {6};
Physical Surface (2) = {10};
// Free surface
Physical Line(1) = {2};
// Bottom of the box.
Physical Line(2) = {3,7};
// Side walls.
Physical Line(3) = {4,9};
// Top (Southmost part)
Physical Line(4) = {1};
