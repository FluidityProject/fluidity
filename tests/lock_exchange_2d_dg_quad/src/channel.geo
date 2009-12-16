Point (1) = {0, 0, 0, 0.01};
Point (2) = {0, 0.1, 0, 0.01};
Line (1) = {1, 2};

// Extrude first half
Extrude {0.4,0,0} {
  Line{1};Layers{40};Recombine;
}
// Extrude second half
Extrude {0.4,0,0} {
  Line{2};Layers{40};Recombine;
}

// Volume number for left of domain.
Physical Surface (1) = {5};
// Volume number for right of domain.
Physical Surface (2) = {9};

// Top of the box.
Physical Line(333) = {4,8};
// Bottom of the box.
Physical Line(444) = {3,7};
// Side walls.
Physical Line(666) = {1,6};
