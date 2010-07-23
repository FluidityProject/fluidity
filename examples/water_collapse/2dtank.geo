Point(1) = {0,0,0,0.1};
Point(2) = {0,0,0,0.1};
Extrude {3.12,0,0} {
  Point{1}; Layers{138};
}
Extrude {0.1,0,0} {
  Point{3}; Layers{30};
}
Extrude {0,0.1,0} {
  Line{1}; Layers{30};
}
Extrude {0,0.1,0} {
  Line{2}; Layers{30};
}
Extrude {0,1.9,0} {
  Line{3}; Layers{82};
}
Extrude {0,1.9,0} {
  Line{7}; Layers{82};
}
Physical Line(19) = {1,2};
Physical Line(20) = {9,17};
Physical Line(21) = {15,11};
Physical Line(22) = {12,4};
Physical Surface(23) = {14,6,10,18};
