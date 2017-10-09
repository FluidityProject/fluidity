Point(1) = {0,0,0,0.1};
Point(2) = {0,1,0,0.1};
Point(3) = {1,1,0,0.1};
Point(4) = {1,0,0,0.1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};

Plane Surface(1) = 1;

Physical Line("left") = 1;
Physical Line("top") = 2;
Physical Line("right") = 3;
Physical Line("bottom") = 4;

Physical Surface("square") = 1;
