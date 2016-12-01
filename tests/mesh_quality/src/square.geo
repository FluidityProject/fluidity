Point (1) = {-1, -1, 0, 0.9};
Point (2) = {-1, 1, 0, 0.9};
Line (1) = {1, 2};

V[]=Extrude {2,0,0} {
  Line{1};
};

Physical Surface(1) = V[1];

Physical Line(1) = {1, V[0], V[2], V[3]};
