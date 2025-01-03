clear all; clc; clf;

load("HEAT.mat");

% hh_x = 1.0 / n_el_x;
% hh_y = 1.0 / n_el_y;
% 
% n_np_x = n_el_x + 1;
% n_np_y = n_el_y + 1;

% [X, Y] = meshgrid(0 : hh_x : 1, 0 : hh_y : 1);
% Z = reshape(disp, n_np_x, n_np_y)';

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
title("Tria numercial solution")

az = -61;
el = 20;
view(az,el);

% EOF