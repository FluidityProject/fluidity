// points and lines:
// p4                   l6                    p3
// --------------------------------------------- 
// |l1  l2                                     |
// --------- p6                                |l5
// p5    l3|                                   |   
//         -------------------------------------
//        p1               l4                 p2 

Point(1) = {5.0,  0.0, 0.0, 0.02};
Point(2) = {25.0 ,0.0 ,0.0, 0.05};
Point(3) = {25.0, 2.0, 0.0, 0.1};
Point(4) = {0.0,  2.0, 0.0, 0.05};
Point(5) = {0.0,  1.0, 0.0, 0.05};
Point(6) = {5.0,  1.0, 0.0, 0.02};


Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};


Line Loop(7) = {1,2,3,4,5,6};
Plane Surface(8) = {7};

// inflow  l1
Physical Line(9) = {4};
// bottom of domain  l2
Physical Line(10) = {5,6,1};
// outflow  l5
Physical Line(11) = {2};
// top  l6
Physical Line(12) = {3};
// Whole domain id
Physical Surface(15) = {8};
