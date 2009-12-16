Point(1) = {-3000000,-10,-4000,10000};
Point(2) = {0,-10,-4000,10000};
Point(3) = {3000000,-10,-4000,10000};
Line(1) = {1,2};
Line(2) = {2,3};
Extrude {0,20,0} {
  Line{1,2}; Layers{1};
}
Extrude {0,0,4000} {
  Surface{6,10}; Layers{1};
}
Physical Volume(55) = {1,2};
//East
Physical Surface(56) = {31};
//West
Physical Surface(57) = {45};
//Top
Physical Surface(58) = {32,54};
//Bottom
Physical Surface(59) = {6,10};
//Front
Physical Surface(60) = {19,41};
//Back
Physical Surface(61) = {27,49};
