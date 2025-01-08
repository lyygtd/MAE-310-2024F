R = 0.5;
L = 1.0;

Point(1) = {L, -L, 0};
Point(2) = {L, L, 0};
Point(3) = {-L, L, 0};
Point(4) = {-L, -L, 0};


Line(1) = {4, 3};
Line(2) = {3, 2};
Line(3) = {2, 1};
Line(4) = {1, 4};


Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};



Transfinite Line{1, 2, 3, 4} = 4;

Transfinite Surface{1};
Transfinite Surface{2};

Recombine Surface{1};


Mesh.ElementOrder = 1;
Mesh.Algorithm = 8;

// EOF

