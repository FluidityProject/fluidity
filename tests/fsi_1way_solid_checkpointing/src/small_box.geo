//

lm = 0.01;

// box points:
Point (6) = {1.0, 0.1, 0, lm};
Point (7) = {1.2, 0.1, 0, lm};
Point (8) = {1.2, 0.2, 0, lm};
Point (9) = {1.0, 0.2, 0, lm};
Line (5) = {9, 6};
Line (6) = {9, 8};
Line (7) = {6, 7};
Line (8) = {8, 7};
// Box surface:
Line Loop(9) = {7, -8, -6, 5};
Plane Surface(10) = {9};




