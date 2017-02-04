size = 0.075;
Point(1) = {0.0,0.0,0,size};
Point(2) = {1.0,0.0,0,size};
Point(3) = {1.0,1.0,0,size};
Point(4) = {0.0,1.0,0,size};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Point(5) = {0.0,0.0,1.0,size};
Point(6) = {1.0,0.0,1.0,size};
Point(7) = {1.0,1.0,1.0,size};
Point(8) = {0.0,1.0,1.0,size};
Line(5) = {8,5};
Line(6) = {5,6};
Line(7) = {6,7};
Line(8) = {7,8};
Line(9) = {8,4};
Line(10) = {5,1};
Line(11) = {2,6};
Line(12) = {3,7};
Line Loop(13) = {5,10,-4,-9};
Plane Surface(14) = {13};
Line Loop(15) = {11,7,-12,-2};
Plane Surface(16) = {15};
Line Loop(17) = {1,11,-6,10};
Plane Surface(18) = {17};
Line Loop(19) = {8,9,-3,12};
Plane Surface(20) = {19};
Line Loop(21) = {5,6,7,8};
Plane Surface(22) = {21};
Line Loop(23) = {1,2,3,4};
Plane Surface(24) = {23};
//left
Physical Surface(25) = {14};
//right
Physical Surface(26) = {16};
//front
Physical Surface(27) = {18};
//back
Physical Surface(28) = {20};
//top
Physical Surface(29) = {22};
//bottom
Physical Surface(30) = {24};
Surface Loop(31) = {22,14,18,24,16,20};
Volume(32) = {31};
Physical Volume(33) = {32};
