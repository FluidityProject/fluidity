
Point (1) = {0, 0, 0, 0};
Extrude {1, 0, 0} {
Point{1};Layers{30};
}
Extrude {0, 1/10, 0} {
Line{1};Layers{10};
}


Physical Line(6) = {3};//left
Physical Line(7) = {4};//right
Physical Line(8) = {2};//up
Physical Line(9) = {1};//down
Physical Surface(10) = {5};
