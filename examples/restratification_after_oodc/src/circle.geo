// Gmsh project created on Fri Feb 08 14:43:37 2008
//radius
// if you want to use fewer resources, make the resolution bigger
res = 5000.0;
Point(1) = {0.0,0.0,0,res};
Point(3) = {250000.0,0.0,0.0,res};
Line(1) = {1,3};


Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{1};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{2};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{5};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{8};
}

Physical Line(14) = {6, 3, 12, 9};//sides
Physical Surface(15) = {4, 7, 10, 13};//surface
