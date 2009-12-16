Point(1) = {0.0,0.0,0.0,600.0};
Point(2) = {32000.0,0.0,0.0,600.0};
Point(3) = {32000.0,32000.0,0.0,600.0};
Point(4) = {0.0,32000.0,0.0,600.0};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {3,4,1,2};
Plane Surface(6) = {5};
Extrude {0,0,-15} {
  Surface{6}; Layers{1};
}
Extrude {0,0,-30} {
  Surface{28}; Layers{1};
}
Extrude {0,0,-60} {
  Surface{50}; Layers{1};
}
Extrude {0,0,-120} {
  Surface{72}; Layers{1};
}
Extrude {0,0,-240} {
  Surface{94}; Layers{1};
}
Extrude {0,0,-480} {
  Surface{116}; Layers{1};
}
Extrude {0,0,-1055} {
  Surface{138}; Layers{1};
}

//top
Physical Surface(161) = {6};
//bottom
Physical Surface(162) = {160};
//east
Physical Surface(163) = {151,129,107,85,63,41,19};
//west
Physical Surface(164) = {159,137,115,93,71,49,27};
//south
Physical Surface(165) = {155,133,111,89,67,45,23};
//north
Physical Surface(166) = {147,125,103,81,59,37,15};

Physical Volume(167) = {7,6,5,4,3,2,1};
