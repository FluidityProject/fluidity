Point (1) = {0, -0.5, 0, 0.025};
Extrude {0.2,0,0} {
  Point{1}; Layers{1}; Recombine;
}
Extrude {0,1.0,0} {
  Line{1}; Layers{20}; Recombine;
}
Line Loop(6) = {3,2,-4,-1};
Physical Line(30) = {4};
Physical Line(32) = {3};
Physical Line(31) = {2};
Physical Line(29) = {1};
Physical Surface(35) = {5};
