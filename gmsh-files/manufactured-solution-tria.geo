R = 0.5;
L = 1.0;

Point(1) = {L, 0, 0};
Point(2) = {L, L, 0};
Point(3) = {0, L, 0};
Point(4) = {0, -L, 0};


Line(1) = {4, 3};
Line(2) = {3, 2};
Line(3) = {2, 1};
Line(4) = {1, 4};


Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};



Transfinite Line{1, 2, 3, 4} = 6;

Transfinite Surface{1};





Mesh.ElementOrder = 1;
Mesh.Algorithm = 8;

// EOF

