Point(1) = {0.4,0.2,0,0.01};
Point(2) = {0.25,0.2,0,0.01};
Point(3) = {-0.3,0.2,0,0.01};
Point(4) = {0.55,0.2,0,0.01};
Point(5) = {1.1,0.2,0,0.01};

Circle(1) = {4,1,2};
Circle(2) = {5,1,3};
Line(3) = {3,2};
Line(4) = {4,5};
Line Loop(5) = {2,3,-1,4};
Plane Surface(6) = {5};
Physical Line(7) = {3};
Physical Line(8) = {1};
Physical Line(9) = {4};
Physical Line(10) = {2};
Physical Surface(11) = {6};
