Point(1) = {0,0,0,0.1};
Extrude {0.4,0,0} {
  Point{1}; Layers{8};
}
Extrude {0,0.4,0} {
  Line{1}; Layers{8};
}
Extrude {0,0.2,0} {
  Line{2}; Layers{4};
}
Extrude {0,0.4,0} {
  Line{6}; Layers{8};
}
Extrude {0.2,0.0,0} {
  Line{12,8,4}; Layers{4};
}
Extrude {0.4,0.0,0} {
  Line{22,18,14}; Layers{8};
}
Physical Line(1) = {1,3,7,11,10,16,36,34,26,27};  // No-slip
Physical Line(2) = {23};                          // Inlet
Physical Line(3) = {30};                          // Outlet
Physical Surface(1) = {13,17,37,9,21,5,29};
Physical Surface(2) = {25};  // Inlet region
Physical Surface(3) = {33};  // Outlet region
