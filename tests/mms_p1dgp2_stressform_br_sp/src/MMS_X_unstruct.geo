pi = 3.141592653589793;
el = pi/XX;

Point(1) = {0.0, 0.0, 0.0, el};
Point(2) = {pi, 0.0, 0.0, el};
Point(3) = {0.0, pi, 0.0, el};
Point(4) = {pi, pi, 0.0, el};
Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 3};
Line(4) = {3, 1};
Line Loop(5) = {2, 3, 4, 1};
Plane Surface(6) = {5};
Physical Line(1) = {4, 3, 2, 1};
Physical Surface(1) = {6};
