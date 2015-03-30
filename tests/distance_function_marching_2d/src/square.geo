Point(1)={0,0,0};
Point(2)={0,1,0};
Point(3)={1,1,0};
Point(4)={1,0,0};

Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};

Line Loop(1)={1,2,3,4};

Plane Surface(1)={1};

Physical Surface(1)={1};


Physical Line(1)={1}; // Left
Physical Line(2)={2}; // Top
Physical Line(3)={3}; // Right
Physical Line(4)={4}; // Bottom

