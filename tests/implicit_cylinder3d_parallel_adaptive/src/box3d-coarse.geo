// box
l1 = 0.03;

Point(1) = {-0.3, -0.3, 0, l1};
Point(2) = {0.3, -0.3, 0, l1};
Point(3) = {-0.3, 0.3, 0, l1};
Point(4) = {0.3, 0.3, 0, l1};

Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 3};
Line(4) = {3, 1};
Line Loop(5) = {2, 3, 4, 1};
Plane Surface(6) = {5};

Extrude {0, 0, -0.1} {
  Surface{6};
Layers{2};
}
Physical Surface(1) = {6}; // top
Physical Surface(2) = {28}; // bottom
Physical Surface(3) = {23}; // inlet
Physical Surface(4) = {15}; // outlet
Physical Surface(5) = {19}; // side
Physical Surface(6) = {27}; // side
Physical Volume(35) = {1};
