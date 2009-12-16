d=-200;
r=10000;

Point(1) = {0,0,0,r};
Point(2) = {200000,0,0,r};
Point(3) = {200000,200000,0,r};
Point(4) = {0,200000,0,r};

Line(1) = {4,1};
Line(2) = {1,2};
Line(3) = {2,3};
Line(4) = {3,4};
Line Loop(5) = {4,1,2,3};
Plane Surface(6) = {5};
Extrude {0,0,d} {
  Surface{6}; Layers{10};
}


Physical Surface(29) = {6};  //top
Physical Surface(30) = {28}; //bottom
Physical Surface(31) = {19}; //east
Physical Surface(32) = {27}; //west
Physical Surface(33) = {23}; //south
Physical Surface(34) = {15}; //north
Physical Volume(35) = {1};
