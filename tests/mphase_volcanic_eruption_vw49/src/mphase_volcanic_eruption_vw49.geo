dx = 1000;
Point(1) = {0, 0, 0, dx};
Point(2) = {150000, 0, 0, dx};
Point(3) = {150400, 0, 0, dx};
Point(4) = {300000, 0, 0, dx};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};

Extrude {0,300000,0} {
Line{1}; Layers{300000/dx};
}
Extrude {0,300000,0} {
Line{2}; Layers{300000/dx};
}
Extrude {0,300000,0} {
Line{3}; Layers{300000/dx};
}

Physical Line(333) = {4, 8, 12};
Physical Line(111) = {5};
Physical Line(222) = {14};
Physical Line(444) = {1, 3};
Physical Line(999) = {2};
Physical Surface(667) = {7,11,15};
