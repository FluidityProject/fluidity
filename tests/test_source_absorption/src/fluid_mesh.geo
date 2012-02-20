// Mesh for the static analysis

// Length of the box
L = 0.4;
// Lenghts of elements
lf = 0.02;
lc = L/(L/lf+1);
// Center line (with different element lengths)
Point(1) = {-0.2,10.2,0,lc};
Point(2) = { 0.2,10.2,0,lc};
// Center line:
Line(1) = {1,2};

// Extruding in Y:
Extrude {0,-0.4,0} {
  Line{1}; Layers{L/lf};
}

// Extruding in Z:
Extrude {0,0,0.2} {
  Surface{5}; Layers{L/(2*lf)};
}
Extrude {0,0,-0.2} {
  Surface{5}; Layers{L/(2*lf)};
}

// Top
Physical Surface(1) = {36, 14};
// Bottom
Physical Surface(2) = {44, 22};
// Sides x
Physical Surface(3) = {48, 26, 40, 18};
// Sides z
Physical Surface(4) = {49, 27};

Physical Volume(5) = {2, 1};
