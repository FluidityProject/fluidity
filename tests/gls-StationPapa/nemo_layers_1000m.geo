// Gmsh project created by create_nemo_column.py
Point(1) = {0,0,0};
Extrude {100.0,0,0} {
  Point{1}; Layers{1};
}
Extrude {0,100.0,0.0} {
  Line{1}; Layers{1}; 
}
//Top
Physical Surface(1) = {5};

col_height = 1000.0;
Extrude {0,0,-col_height} {
	Surface{5}; Layers{ {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}, {0.003,0.0091,0.0154,0.0219,0.0286,0.0355,0.0427,0.0502,0.0581,0.0664,0.0753,0.0847,0.0948,0.1057,0.1176,0.1305,0.1446,0.1601,0.1773,0.1963,0.2174,0.2409,0.2671,0.2964,0.3292,0.3659,0.4068,0.4525,0.5034,0.56,0.6227,0.6919,0.768,0.8514,0.9422,1.00}};
}
Physical Line(44) = {3};
Physical Line(45) = {2};
Physical Line(46) = {4};
Physical Line(47) = {1};
Physical Volume(100) = {1};
Physical Surface(101) = {18};
Physical Surface(102) = {22};
Physical Surface(103) = {26};
Physical Surface(104) = {14};
Physical Surface(2) = {27};
