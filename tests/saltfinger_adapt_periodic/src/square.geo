Point(1) = {0,0,0,0.1};
Extrude {1,0,0} {
  Point{1}; Layers{100};Recombine;
}
Extrude {0,2,0} {
  Line{1}; Layers{100};Recombine;
}
Physical Line(6) = {2};//t
Physical Line(7) = {1};//b
Physical Line(8) = {3};//l
Physical Line(9) = {4};//r
Physical Surface(10) = {5};
