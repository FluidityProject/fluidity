Point(1) = {0,0,0,0.1};
Extrude {1,0,0} {
  Point{1};Layers{1};
}
Extrude {0,1,0} {
  Line{1};Layers{1};
}
Extrude {0,0,1} {
  Surface{5};Layers{1};
}
