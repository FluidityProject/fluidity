Point(1) = {-10,0,0,0.01};
Point(2) = {10,0,0,0.01};
Point(3) = {-10.0,1.0,0,0.01};
Point(4) = {10.0,1.0,0,0.01};
Line(1) = {1,2};
Line(2) = {2,4};
Line(4) = {3,1};
Line(5) = {4,3};
Line Loop(8) = {1,2,4,5};

Plane Surface(9) = {8};
Physical Surface(10) = {9};
