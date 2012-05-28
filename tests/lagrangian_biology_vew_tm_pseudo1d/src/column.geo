Point(1) = {0,0,0};
Extrude {10.0,0,0} {
  Point{1}; Layers{1};
}
Extrude {0,10.0,0.0} {
  Line{1}; Layers{1}; 
}
Extrude {0.0,0.0,-300.0} {
  Surface{5}; Layers{300}; 
}
//Top
Physical Surface(38) = {5};
//Bottom
Physical Surface(39) = {27};
//North
Physical Surface(40) = {22};
//South
Physical Surface(41) = {14};
//West
Physical Surface(42) = {26};
//East
Physical Surface(43) = {18};
Physical Volume(44) = {1};

