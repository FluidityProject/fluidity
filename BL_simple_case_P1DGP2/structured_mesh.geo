c = 3;

Point(1) = {0., 0, 0.};
Point(2) = {1, 0, 0.};
Point(3) = {1, 1/c, 0.};
Point(4) = {0., 1/c, 0.};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(5) = {3, 4, 1, 2};
Ruled Surface(6) = {5};

Transfinite Line {1, 3} = c*2+1 Using Progression 1;
Transfinite Line {-2, 4} = 3 Using Progression 1;

Transfinite Surface {6} Alternate;

Physical Line(7) = {4};//left
Physical Line(8) = {2};//right
Physical Line(9) = {3};//up
Physical Line(10) = {1};//down
Physical Surface(11) = {6};
