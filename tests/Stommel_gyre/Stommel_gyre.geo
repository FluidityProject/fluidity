Point (1) = {0.0, 0.0, 0.0, 0.02};
Point (2) = {0.0, 1.0, 0.0, 0.02};
Line (1) = {1, 2};

//Extrude to make box
Extrude {1.0,0,0} {
  Line{1}; Layers{50};
}
// top
Physical Line(6) = {4};
// bottom
Physical Line(7) = {3};
// ends
Physical Line(8) = {2,1};
Physical Surface(9) = {5};
