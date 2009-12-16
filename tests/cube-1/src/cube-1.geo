Point(1) = {0,0,0,1.0};
Point(2) = {1,0,0,1.0};
Line(1) = {1,2};
Extrude {0,1,0} {
  Line{1};
}
Extrude {0,0,1} {
  Surface{5};
}
Physical Volume(1) = {1};
//Sides
Physical Surface(1) = {14,22};
//TopBottom
Physical Surface(2) = {5,27};
//Inflow
Physical Surface(3) = {26};
//Outflow
Physical Surface(4) = {18};
