Point(1) = {0,0,0,1e5};
Extrude {1e6,0,0} {
  Point{1};
}
Extrude {0,2e5,0} {
  Line{1};
}
// Reserve 1 and 2 for top and bottom of extruded mesh
// Outer ends of the channel (x=0 and x=1e6)
Physical Line(3) = {3};
Physical Line(4) = {4};
// Sides of the channel (y=0 and y=1e5)
Physical Line(5) = {1};
Physical Line(6) = {2};
// 
Physical Surface(1) = {5};
