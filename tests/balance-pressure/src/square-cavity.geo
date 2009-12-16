Point (1) = {0, 0, 0, 0.025};
Point (2) = {0, 1, 0, 0.025};
Line (1) = {1, 2};

Extrude {1,0,0} {
  Line{1};
}

// Face
Physical Surface (1) = {5};
// Top
Physical Line(2) = {4};
// Bottom
Physical Line(3) = {3};
// Left
Physical Line(4) = {1};
// Right
Physical Line(5) = {2};
