// Geometry parameters [m]
length=0.2;
height=0.2;

Point(1) = {0, 0, 0, 1.0};
Point(2) = {length, 0, 0, 1.0};
Point(3) = {length, height, 0, 1.0};
Point(4) = {0, height, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {3, 4, 1, 2};

Plane Surface(1) = {1};
//+
Physical Curve("bottom_edge") = {1};
//+
Physical Curve("right_edge") = {2};
//+
Physical Curve("top_edge") = {3};
//+
Physical Curve("left_edge") = {4};
//+
Physical Surface("body") = {1};
