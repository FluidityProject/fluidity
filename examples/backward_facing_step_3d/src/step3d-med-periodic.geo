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

Point(1) = {30, 0, 1, 0.25};
//Inlet
Extrude {0, 0, 5} {
  Point{1}; Layers{20};
}
Extrude {0, 0, -1} {
  Point{1}; Layers{4};
}

//Extrude along
Extrude {-40, 0, 0} {
  Line{1}; Layers{160};
}
Extrude {-40, 0, 0} {
  Line{2}; Layers{160};
}
Delete {
  Surface{6,10};
  Line{7,4,5,9};
  Point{7};
}

//Draw new lines
Point(8) = {0, 0, 0, 0.25};
Extrude {0, 0, 1} {
  Point{8}; Layers{4};
}
Extrude {30, 0, 0} {
  Point{8}; Layers{120};
}
Extrude {10, 0, 0} {
  Point{4}; Layers{40};
}
Line(13) = {5,2};
Line Loop(14) = {12,-10,11,-2,1,-13,-3};

Plane Surface(5) = {14};

//Extrude volume
Extrude {0, 4, 0} {
  Surface{5}; Layers{16};
}

// in
Physical Surface(52) = {50};
//out
Physical Surface(53) = {38,42};
//top
Physical Surface(54) = {46};
//front
Physical Surface(55) = {5};
//back
Physical Surface(56) = {51};
//bottom of domain
Physical Surface(57) = {26,30,34};

Physical Volume(58) = {1};

