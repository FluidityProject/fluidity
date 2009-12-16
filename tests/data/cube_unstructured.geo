Point(1) = {0,0,0,0.1};
Point(2) = {1,0,0,0.1};
Line(1) = {1,2};
Extrude {0,1,0} {
  Line{1};
}
Extrude {0,0,1} {
  Surface{5};
}
Physical Volume(1) = {1};
//Top
Physical Surface(1) = {27};
//Bottom
Physical Surface(2) = {5};
//Sides
Physical Surface(3) = {14,22,26,18};
