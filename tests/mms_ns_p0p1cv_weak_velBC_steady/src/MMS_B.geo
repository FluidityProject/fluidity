x = 2.22144146908;
c = 1.57079632679;

Point(1) = {x+c,0.0+c,0,0.32};
Point(2) = {0.0+c,x+c,0,0.32};
Point(3) = {-x+c,0.0+c,0,0.32};
Point(4) = {0.0+c,-x+c,0,0.32};
Line(1) = {4,3};
Line(2) = {3,2};
Line(3) = {2,1};
Line(4) = {1,4};
Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};
Physical Line(7) = {3};
Physical Line(8) = {2};
Physical Line(9) = {1};
Physical Line(10) = {4};
Physical Surface(12) = {6};
