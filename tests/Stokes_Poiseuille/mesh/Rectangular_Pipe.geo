// Define rectangular half-height
h = 1;
// Define rectangular length
l = 2e1;
// Define the mesh resolution
r = 4e-1;
// Define number of points on the horizontal lines
nh = Round(l / r);
// Define number of points on the vertical lines
nv = Round(2 * h / r);

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

Transfinite Line{1, 3} = nh;
Transfinite Line{2, 4} = nv;
Transfinite Surface{1} = {1, 2, 3, 4};
