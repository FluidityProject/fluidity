nlayers = 20;
Point (1) = {0, 0, 0, 1.0/nlayers};
Point (2) = {1, 0, 0, 1.0/nlayers};
Line (1) = {1, 2};

// Extrude first half
Extrude {0,1,0} {
  Line{1};Layers{nlayers};Recombine;
}
Physical Line(6) = {3};
Physical Line(7) = {1};
Physical Line(8) = {4};
Physical Line(9) = {2};
Physical Surface(10) = {5};
