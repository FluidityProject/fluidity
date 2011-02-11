Point(1) = {0, 0, 0, 1e+22};
Point(2) = {360, 0, 0, 1e+22};
Line(1) = {1, 2};

Extrude {0,360,0} {
Line{1}; Layers{10}; Recombine;
}