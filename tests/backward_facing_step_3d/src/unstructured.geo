// points and lines
//
//
//        ----------------------------------------------
//       / |      |                                    /|
//      /                                             / |
//     /   |      |                                  /  |
//    /       --                                    /   |
//  p4----/------/----------------------------------p3  |
//    |           |                                |    |
//    | /      /   --  --  --  --  --  --  --   -- | - / 
//    |         /                                  |  /
//    --------p6                                   | /
//    p5     |/                                    |/   
//           ---------------------------------------
//           p1                                    p2        


// Values based on Le, Moin and Kim, 1997, JFM 330, 349-374.
Point(1) = {0.0,   0.0, 0.0, 0.1};
Point(2) = {30.0,  0.0 ,0.0, 0.2};
Point(3) = {30.0,  0.0, 6.0, 0.2};
Point(4) = {-10.0, 0.0, 6.0, 0.05};
Point(5) = {-10.0, 0.0, 1.0, 0.05};
Point(6) = {0.0,   0.0, 1.0, 0.05};


Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};


Line Loop(7) = {5,6,1,2,3,4};
Plane Surface(8) = {7};
Extrude {0,4,0} {
  Surface{8};
}

// in
Physical Surface(41) = {39};
//out
Physical Surface(42) = {31};
//step_down
Physical Surface(43) = {23};
//bottom_step
Physical Surface(44) = {19};
//bottom_main
Physical Surface(45) = {27};
//top
Physical Surface(46) = {35};
//front
Physical Surface(47) = {8};
//back
Physical Surface(48) = {40};

Physical Volume(49) = {1};
