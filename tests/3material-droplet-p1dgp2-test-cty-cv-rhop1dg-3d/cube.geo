// define a layer variable
lay_x = 5;
lay_y = 5;
lay_z = 20;

// define a len variable
len_x = 1.0;
len_y = 1.0;
len_z = 4.0;

Point(1) = {0.0, 0.0, 0.0, 1.0};

Extrude {len_x, 0.0, 0.0} {
  Point{1}; Layers{lay_x}; 
}

Extrude {0.0, len_y, 0.0} {
  Line{1}; Layers{lay_y}; 
}

Extrude {0.0, 0.0, len_z} {
  Surface{5}; Layers{lay_z}; 
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
