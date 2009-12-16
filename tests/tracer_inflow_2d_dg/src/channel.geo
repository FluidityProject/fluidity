Point (1) = {0, -1, 0, 1.9};
Point (2) = {0, 1, 0, 1.9};
Line (1) = {1, 2};

Extrude {10,0,0} {
  Line{1};Layers{20};
}

// Volume number for whole domain.
Physical Surface (1) = {5};
// Top and bottom of the box.
Physical Line(333) = {3,4};
// Right wall
Physical Line(444) = {2};
// Left wall
Physical Line(666) = {1};
