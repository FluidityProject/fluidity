Point(1) = {-24.0,-8.0,0.0,0.4};
Point(2) = {-24.0,8.0,0.0,0.4};
Point(3) = {24.0,8.0,0.0,0.4};
Point(4) = {24.0,-8.0,0.0,0.4};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};


Line Loop(5) = {2, 3, 4, 1};
Plane Surface(6) = {5};
//north
Physical Line(7) = {2};
//south
Physical Line(8) = {4};
//west
Physical Line(9) = {1};
//east
Physical Line(10) = {3};
Physical Surface(11) = {6};
