dx = 100;
Point(1) = {0, 0, 0, dx};
Point(2) = {200, 0, 0, dx};
Point(3) = {7000, 0, 0, dx};
Line(1) = {1, 2};
Line(2) = {2, 3};

Extrude {0,7000,0} {
Line{1}; Layers{7000/dx};
}
Extrude {0,7000,0} {
Line{2}; Layers{7000/dx};
}

Physical Line(333) = {3, 7};
Physical Line(111) = {4};
Physical Line(222) = {9};
Physical Line(444) = {2};
Physical Line(999) = {1};
Physical Surface(667) = {6,10};
