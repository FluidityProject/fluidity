Point (1) = {0.0, 0, 0, 0.001};
Point (2) = {0.0, 0.001, 0, 0.001};
Line (1) = {1, 2};

//Extrude to make box
Extrude {0.8,0,0} {
  Line{1}; Layers{80};
}
// Extrude up
Extrude {0,0,0.1} {
  Surface{9,5}; Layers{10};
}
// top
Physical Surface(28) = {27};
// bottom
Physical Surface(29) = {5};
// E-W (x)
Physical Surface(30) = {22,14};
// N-S (y)
Physical Surface(31) = {18,26};
Physical Volume(32) = {1};
