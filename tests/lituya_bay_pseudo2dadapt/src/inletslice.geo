Point(1) = {122,-122,0,100};
Point(2) = {1220,-122,0,100};
Point(3) = {2000,658,0,100};
Point(4) = {2000,1100,0,100};
Point(5) = {-1100,1100,0,100};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,1};
Line Loop(6) = {1,2,3,4,5};
Plane Surface(7) = {6};
Extrude {0,0,1} {
  Surface{7}; Layers{2};
}
Physical Surface(8) = {17};
Physical Surface(9) = {21};
Physical Surface(10) = {25};
Physical Surface(11) = {29};
Physical Surface(12) = {33};
Physical Surface(13) = {34};
Physical Surface(14) = {7};
Physical Volume(42) = {1};
