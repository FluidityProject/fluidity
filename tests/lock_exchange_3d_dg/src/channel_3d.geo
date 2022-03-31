Point (1) = {0.0, 0, 0, 0.001};
Point (2) = {0.0, 0.001, 0, 0.001};
Line (1) = {1, 2};

//Extrude to make box
Extrude {0.4,0,0} {
  Line{1}; Layers{30};
}
Extrude {0.4,0,0} {
  Line{2}; Layers{30};
}

Extrude {0,0,0.1} {
  Surface{5,9}; Layers{10};
}
// Sides
Physical Surface(111) = {48,18};
// Top and bottom
Physical Surface(222) = {53,5,31,9};
// Front and back
Physical Surface(333) = {44,52,22,30};
// Left
Physical Volume(1) = {1};
// Right
Physical Volume(2) = {2};
