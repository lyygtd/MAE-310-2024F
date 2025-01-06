clear all; clc; clf;

load("Solid-quad.mat");

pos = 4; 

logicalIdx = true(size(x_coor));
logicalIdx(pos) = false;
 

x = x_coor(logicalIdx);
y = y_coor(logicalIdx);
z = disp(logicalIdx);


scatter3(x, y, z, 6, z, 'filled');

xlabel("x")
ylabel("y")
zlabel("z")
shading interp
title("Quad numercial solution")

clear

load("Solid-tria.mat");

pos = 4; 

logicalIdx = true(size(x_coor));
logicalIdx(pos) = false;


x = x_coor(logicalIdx);
y = y_coor(logicalIdx);
z = disp(logicalIdx);

figure
scatter3(x, y, z, 6, z, 'filled');

xlabel("x")
ylabel("y")
zlabel("z")
shading interp
title("Tria numercial solution")
% EOF