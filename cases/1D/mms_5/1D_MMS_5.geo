// Gmsh project for mms_5 (same geometry as mms_4)
SetFactory("OpenCASCADE");

Point(1) = {-4, 0, 0, 1.0};
Point(2) = {0, 0, 0, 1.0};
Point(3) = {4, 0, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};

Physical Point("bc:west", 1) = {1};
Physical Point("bc:east", 2) = {3};

Physical Curve("material:mms_material_1", 4) = {1};
Physical Curve("material:mms_material_2", 5) = {2};

Physical Curve("source:mms_source_1",6) = {1};
Physical Curve("source:mms_source_2", 7) = {2};

