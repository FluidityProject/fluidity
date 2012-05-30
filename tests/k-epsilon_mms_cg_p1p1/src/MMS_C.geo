pi = 3.1415926535897931;
e = 40;

Point(1) = {0.0, 0.0, 0.0, 1};

Extrude {pi, 0.0, 0.0} {
  Point{1}; Layers{e};
}

Extrude {0.0, pi, 0.0} {
  Line{1}; Layers{e};
}

Physical Line(7) = {1};  //reentrainment
Physical Line(9) = {2};  //lid
Physical Line(10) = {3};  //in
Physical Line(8) = {4};  //out
Physical Surface(12) = {5};
