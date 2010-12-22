Point(1) = {0, 0, 0};
Point(2) = {5.488, 0, 0};
Point(3) = {5.488, 3.402, 0};
Point(4) = {0, 3.402, 0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {3, 4, 1, 2};
Plane Surface(6) = {5};
Physical Line(7) = {3};
Physical Line(8) = {1};
Physical Line(9) = {2};
Physical Line(10) = {4};
Physical Surface(11) = {6};

Field[1] = MathEval;
//Field[1].F = "100" // Coarse
//Field[1].F = "((x-4.521)*(x-4.521)+(y-1.696)*(y-1.696)^0.5)/10+0.1"; // Quadratic around gauge 2
//Field[1].F = "((x-4.521)*(x-4.521))^0.5/10+0.1"; // linear refining towards right
//Field[1].F = "((x-4.521)*(x-4.521))^0.5/10+0.1"; // linear refining towards right plus refinment at the island
Field[1].F = "min(((x-4.521)*(x-4.521))^0.5/10 + 0.1, (((x-3.39)^2+(y-1.659)^2)^0.5)^4/1+0.05)"; // linear refining towards right plus refinment at the island
Background Field = 1;

