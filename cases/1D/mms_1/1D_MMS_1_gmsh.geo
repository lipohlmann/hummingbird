// Gmsh project created on Thu Sep  3 14:05:13 2026
SetFactory("OpenCASCADE");
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(1) = {-1, 0, 0, 1.0};
//+
Line(1) = {2, 1};
//+
Physical Point("bc:west", 2) = {1};
//+
Physical Point("bc:east", 3) = {2};
//+
Physical Curve("material:mms_material", 4) = {1};
//+
Physical Curve("source:mms_source", 5) = {1};
