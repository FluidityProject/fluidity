// Gmsh project created on Thu Mar  6 16:27:33 2008
Point(1) = {0,0,0,10};
Extrude {19000-(800/(Tan(2.5*Pi/180))),0,0} {
  Point{1}; 
}
Extrude {800/(Tan(2.5*Pi/180)),0,0} {
  Point{2};
}
Extrude {500,0,0} {
  Point{3}; 
}
Extrude {500,0,0} {
  Point{4}; 
}
Extrude {0,1000,0} {
  Line{1,2,3,4}; Layers{100};
}
// top
Physical Line(21) = {17,13,9,5};
//bottom
Physical Line(22) = {4,3,2,1};
// x = 0
Physical Line(23) = {6};
// x = 20km
Physical Line(24) = {19};
// plateau
Physical Surface(25) = {20,16};
// rest of domain
Physical Surface(26) = {12,8};
