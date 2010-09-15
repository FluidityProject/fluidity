Point (1) = {0.0, 0, 0, 0.005};
Point (2) = {0.0, 0.1, 0, 0.005};
Line (1) = {1, 2};

//Extrude to make box,
//putting extra resolution in near interface
Extrude {0.395,0,0} {
  Line{1}; Layers{10};
}
Extrude {0.01,0,0} {
  Line{2}; Layers{50};
}
Extrude {0.395,0,0} {
  Line{6}; Layers{10};
}
// top
Physical Line(6) = {4,8,12};
// bottom
Physical Line(7) = {3,7,11};
// ends
Physical Line(8) = {1,10};
Physical Surface(9) = {5,9,13};
