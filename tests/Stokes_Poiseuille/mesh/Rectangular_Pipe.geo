// Define rectangular half-height
h = 1;
// Define rectangular length
l = 2e1;

// Define the four corners of the rectangle, clockwise starting from the
// top-left one
Point(1) = {0, h, 0};
Point(2) = {l, h, 0};
Point(3) = {l, -h, 0};
Point(4) = {0, -h, 0};

// Define the boundary lines of the rectangular domain, clockwise starting from
// the top-left corner
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
// Define the surface delimited by the lines
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Physical Surface(1) = {1};

// Set the elements sizes : a box is defined inside the domain and the
// resolution of the mesh both inside and outside the box is set to be the same
Field[1] = Box;
Field[1].VIn = 3e-1;
Field[1].VOut = 3e-1;
Field[1].XMax = 3 * l / 4;
Field[1].XMin = l / 4;
Field[1].YMax = h / 2;
Field[1].YMin = -h / 2;

Background Field = 1;

// Force triangular elements
Mesh.RecombineAll = 0;
// Prevent elements sizes to be modified close to the boundaries
Mesh.CharacteristicLengthExtendFromBoundary = 0;
