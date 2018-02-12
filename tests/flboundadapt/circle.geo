Point(1) = {-1,-1,0};
Point(2) = { 1,-1,0};
Point(3) = { 1, 1,0};
Point(4) = {-1, 1,0};

Point(5) = {0,0,0};

r=0.3;

Point(6) = {-r, 0,0};
Point(7) = { 0, r,0};
Point(8) = { r, 0,0};
Point(9) = { 0,-r,0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Circle(5) = {6,5,7};
Circle(6) = {7,5,8};
Circle(7) = {8,5,9};
Circle(8) = {9,5,6};

Line Loop(1) = {1,2,3,4};
Line Loop(2) = {5,6,7,8};

Plane Surface(1) = {1,2};

Physical Surface(1) = {1};

Physical Line(1) = {4};
Physical Line(2) = {1,3};
Physical Line(4) = {2};
Physical Line(5) = {5,6,7,8};

Field[1]=MathEval;
Field[1].F="0.2";

Background Field =2;