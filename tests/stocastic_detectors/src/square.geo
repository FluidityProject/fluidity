r=0.05;
Point (1) = {-0.5, -0.5, 0,r};
Point (2) = {-0.5, 0.5, 0,r};
Line (1) = {1, 2};

Extrude {1,0,0} {
  Line{1};
}

// Volume number for whole domain.
Physical Surface (1) = {5};
// Top of the box.
Physical Line(333) = {4};
// Bottom of the box.
Physical Line(444) = {3};
// Side walls.
Physical Line(666) = {1,2};
