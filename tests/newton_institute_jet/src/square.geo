a = 6e6;
N = 50.;
dx = 2 * a / N;
Point (1) = {0, 0, 0, dx};
Point (2) = {0, 2*a, 0, dx};
Line (1) = {1, 2};

Extrude {2*a,0,0} {
  Line{1};Layers{N};
}

// Volume number for whole domain.
Physical Surface (1) = {5};
// Top of the box.
Physical Line(4) = {4};
// Bottom of the box.
Physical Line(3) = {3};
// Left
Physical Line(1) = {1};
// Right
Physical Line(2) = {2};
