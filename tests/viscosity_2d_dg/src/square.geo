Point (1) = {-1, -1, 0, 0.9};
Point (2) = {-1, 1, 0, 0.9};
Line (1) = {1, 2};

Extrude {1,0,0} {
  Line{1};Layers{1};
}

// Volume number for whole domain.
Physical Surface (1) = {5};
// Top of the box.
Physical Line(333) = {4};
// Bottom of the box.
Physical Line(444) = {3};
// Side walls.
Physical Line(666) = {1,2};
