// define a layer variable
lay = 5;

// define a len variable
len = 300.0;

Point(1) = {0.0, 0.0, 0.0, 1.0};

Extrude {len, 0.0, 0.0} {
  Point{1}; Layers{lay}; 
}

Extrude {0.0, len, 0.0} {
  Line{1}; Layers{lay}; 
}

Extrude {0.0, 0.0, len} {
  Surface{5}; Layers{lay}; 
}

// bottom
Physical Surface(1) = {5};

// top
Physical Surface(2) = {27};

// left
Physical Surface(3) = {26};

// right
Physical Surface(4) = {18};

// front
Physical Surface(5) = {14};

// back
Physical Surface(6) = {22};

Physical Volume(1) = {1};
