Point(1) = {0.1,-0.3,0,0.02};
Point(2) = {0.6,-0.3,0,0.02};
Line (1) = {1, 2};
Extrude {0,0.4,0} {
  Line{1}; Layers{20};
}
// y const
Physical Line(7) = {1};
// x const
Physical Line(10) = {3};
// y const
Physical Line(9) = {2};
// x const
Physical Line(8) = {4};
Physical Surface(11) = {5};
