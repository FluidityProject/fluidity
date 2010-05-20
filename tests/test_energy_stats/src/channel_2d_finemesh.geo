Point (1) = {0.0, 0, 0, 0.0005};
Point (2) = {0.0, 0.1, 0, 0.0005};
Line (1) = {1, 2};

//Extrude to make box
Extrude {0.8,0,0} {
  Line{1}; Layers{1600};
}
// top
Physical Line(6) = {4};
// bottom
Physical Line(7) = {3};
// ends
Physical Line(8) = {2,1};
Physical Surface(9) = {5};
