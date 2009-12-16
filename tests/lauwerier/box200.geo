Point(1) = {0,0,0,0.1};
Point(2) = {3.14,0,0,0.1};
Point(3) = {3.14,6.28,0,0.1};
Point(4) = {0.,6.28,0,0.1};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};
Extrude {0,0,0.00157} {
  Surface{6}; Layers{1};
}
Physical Line(29) = {10,3};
Physical Line(30) = {22,18,9,2,14,1,8,13,4,11};
Physical Surface(31) = {28};
// top
Physical Surface(32) = {6};
// bottom
Physical Surface(33) = {23};
//north
Physical Surface(34) = {27,19};
//east, west (resp)
Physical Surface(35) = {15};
//south
Physical Volume(36) = {1};
