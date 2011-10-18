Point(1) = {0,0,0,100};
Point(2) = {100,0,0,100};
Point(3) = {100,100,0,100};
Point(4) = {0,100,0,100};
Line(1) = {4,3};
Line(2) = {3,2};
Line(3) = {2,1};
Line(4) = {1,4};
Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};
Extrude {0,0,-200} {
  Surface{6};Layers{20};
}
cells[0] = 1;
heights[0] = .0333;
thickness[0]=.0333;
For i In {1:9}
  cells[i] = 1;
  thickness[i]=thickness[i-1]*1.18;
  heights[i] = heights[i-1]+thickness[i];
EndFor
cells[10] = 1;
heights[10] = 1;


Extrude {0,0,-300} {
  Surface{28}; Layers{ cells[], heights[] };
}

Extrude {0,0,-2500} {
  Surface{50};
}
Physical Surface(111) = {6};
Physical Surface(222) = {27,15,19,23,41,45,49,37,59,63,67,71};
Physical Surface(333) = {72};
Physical Volume(1) = {1};
Physical Volume(2) = {2};
Physical Volume(3) = {3};
