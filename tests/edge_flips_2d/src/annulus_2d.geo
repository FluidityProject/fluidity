Point(1) = { 0, 0, 0};

Point(2) = { 1, 0, 0};
Point(3) = { 0, 1, 0};
Point(4) = {-1, 0, 0};
Point(5) = { 0,-1, 0};


Point(6) = { 2, 0, 0};
Point(7) = { 0, 2, 0};
Point(8) = {-2, 0, 0};
Point(9) = { 0,-2, 0};


Point(10) = { 3, 0, 0};
Point(11) = { 0, 3, 0};
Point(12) = {-3, 0, 0};
Point(13) = { 0,-3, 0};

Circle(1)={2,1,3};
Circle(2)={3,1,4};
Circle(3)={4,1,5};
Circle(4)={5,1,2};

Circle(5)={6,1,7};
Circle(6)={7,1,8};
Circle(7)={8,1,9};
Circle(8)={9,1,6};

Circle(9) ={10,1,11};
Circle(10)={11,1,12};
Circle(11)={12,1,13};
Circle(12)={13,1,10};


Line Loop(1)={1,2,3,4};
Line Loop(2)={5,6,7,8};
Line Loop(3)={9,10,11,12};

Plane Surface(1) = {3,2};
Plane Surface(2) = {2,1};

Physical Surface(1)={1};
Physical Surface(2)={2};

Physical Line(2)={1,2,3,4};
Physical Line(1)={9,10,11,12};