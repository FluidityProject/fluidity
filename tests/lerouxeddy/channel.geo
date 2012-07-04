Ls=10;
Lsx=Ls*5;
Lsy=Ls*3;
Lx=2000000;
Ly=1200000;
Point(1) = {-Lx*0.5, -Ly*0.5, 0, 0.1*Lx};
Extrude {Lx, 0, 0} {
  Point{1};Layers{Lsx};
}
Extrude {0., Ly, 0} {
  Line{1}; Layers{Lsy};
}
Physical Line(1) = {3};
Physical Line(2) = {4};
Physical Line(3) = {1};
Physical Line(4) = {2};
Physical Surface(5) = {5};
