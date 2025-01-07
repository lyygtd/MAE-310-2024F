clear all; clc;

kappa = 1.0; % conductivity


% exact solution
exact = @(x,y) x*(1-x)*y*(1-y);
exact_x = @(x,y) (1-2*x)*y*(1-y);
exact_y = @(x,y) x*(1-x)*(1-2*y);

f = @(x,y) 0; % source term

% quadrature rule
n_int = 3;
weight = [1/3, 1/3, 1/3]/2;
xi =     [2/3, 1/6, 1/6];
eta =    [1/6, 2/3, 1/6];

% quadrature rule for h boundary
n_int_h = 10;
[xi_h, weight_h] = Gauss(n_int_h, -1, 1);

% mesh数据导入
fileID = fopen("..\..\gmsh-files\quarter-plate-with-hole-tria.msh", 'r');
if fileID == -1
    error('无法打开文件');
end

x_coor = [];
y_coor = [];
maxCols = 0;
inDataSection = false;

tline = fgetl(fileID);

%生成x_coor和y_coor
while ischar(tline)
    if strcmp(tline, '$Nodes')
        inDataSection = true;
        tline = fgetl(fileID);
        continue;
    elseif strcmp(tline, '$EndNodes')
        inDataSection = false;
    elseif inDataSection
        rowStrs = strsplit(tline);
        rowData = str2double(rowStrs);
        numRows = numel(rowData);
        if numRows > maxCols
            maxCols = numRows;
        end
        rowDataPadded = [rowData; zeros(maxCols - numRows, 1)]';
        rowDataPadded = rowDataPadded(1:maxCols);
        if numRows == 4
            x_coor = [x_coor; rowDataPadded(end-2)'];
            y_coor = [y_coor; rowDataPadded(end-1)'];
        end
    end
    tline = fgetl(fileID);
end


fileID = fopen("..\..\gmsh-files\quarter-plate-with-hole-tria.msh", 'r');

IEN = [];
n_el = 0;
maxCols = 0;
inDataSection = false;

tline = fgetl(fileID);

%生成IEN
while ischar(tline)
    if strcmp(tline, '$Elements')
        inDataSection = true;
        tline = fgetl(fileID);
        continue;
    elseif strcmp(tline, '$EndElements')
        inDataSection = false;
    elseif inDataSection
        rowStrs = strsplit(tline);
        rowData = str2double(rowStrs);
        numRows = numel(rowData);
        if numRows > maxCols
            maxCols = numRows;
        end
        rowDataPadded = [rowData; zeros(maxCols - numRows, 1)]';
        rowDataPadded = rowDataPadded(1:maxCols);
        if numRows == 8
            if rowDataPadded(5) == 1
                temp_ien = rowDataPadded(end-2 : end)';
                sj = temp_ien(3);
                temp_ien(3) = temp_ien(2);
                temp_ien(2) = sj;
                IEN = [IEN; temp_ien];
                n_el = n_el + 1;
            else
                IEN = [IEN; rowDataPadded(end-2 : end)'];
                n_el = n_el + 1;
            end
        end
    end
    tline = fgetl(fileID);
end


fclose(fileID);
n_en   = 3;               % number of nodes in an element
n_np   = length(x_coor); % total number of nodal points


gg = zeros(n_np,1); % g boundary condition
left_h_boundary_nodes = [];
right_h_boundary_nodes = [];
bottom_h_boundary_nodes = [];
top_h_boundary_nodes = [];
circle_h_boundary_nodes = [];

% Boundary condition settings

left_Boudary_is_Dirihelt = true;
g_left = @(y) 0;
left_Boudary_is_Neummen = false;
h_left = @(y) 40;

right_Boudary_is_Dirihelt = false;
g_right = @(y) 0;
right_Boudary_is_Neummen = true;
h_right = @(y) y;

top_Boudary_is_Dirihelt = false;
g_top = @(x) x;
top_Boudary_is_Neummen = true;
h_top = @(x) x;

bottom_Boudary_is_Dirihelt = true;
g_bottom = @(x) 0;
bottom_Boudary_is_Neummen = false;
h_bottom = @(x) 0;

circle_Boudary_is_Dirihelt = false;
g_circle = @(x,y) 0;
circle_Boudary_is_Neummen = true;
h_circle = @(x,y) 0;

% ID array
ID = zeros(n_np,1);
counter = 0;

for nn = 1 : n_np
    if x_coor(nn) == -1   % left boundary
        if left_Boudary_is_Neummen
            counter = counter +1;
            ID(nn) = counter;
            left_h_boundary_nodes = [left_h_boundary_nodes; nn];
        elseif left_Boudary_is_Dirihelt
            ID(nn) = 0;
            gg(nn) = g_left(y_coor(nn));
        end

    elseif x_coor(nn) == 1  % right boundary
        if right_Boudary_is_Neummen
            counter = counter +1;
            ID(nn) = counter;
            right_h_boundary_nodes = [right_h_boundary_nodes; nn];
        elseif right_Boudary_is_Dirihelt
            ID(nn) = 0;
            gg(nn) = g_right(y_coor(nn));
        end

    elseif y_coor(nn) == 1  % top boundary
        if top_Boudary_is_Neummen
            counter = counter +1;
            ID(nn) = counter;
            top_h_boundary_nodes = [top_h_boundary_nodes; nn];
        elseif top_Boudary_is_Dirihelt
            ID(nn) = 0;
            gg(nn) = g_top(x_coor(nn));
        end

    elseif y_coor(nn) == -1  % bottom boundary
        if bottom_Boudary_is_Neummen
            counter = counter +1;
            ID(nn) = counter;
            bottom_h_boundary_nodes = [bottom_h_boundary_nodes; nn];
        elseif bottom_Boudary_is_Dirihelt
            ID(nn) = 0;
            gg(nn) = g_bottom(x_coor(nn));
        end

    elseif abs(sqrt((x_coor(nn)+1)^2 + (y_coor(nn)+1)^2) - 0.5) <= 1e-14 % circle boundary
        if bottom_Boudary_is_Neummen
            counter = counter +1;
            ID(nn) = counter;
            circle_h_boundary_nodes = [circle_h_boundary_nodes; nn];
        elseif bottom_Boudary_is_Dirihelt
            ID(nn) = 0;
            gg(nn) = g_circle(x_coor(nn), y_coor(nn));
        end

    else
        counter = counter +1;
        ID(nn) = counter;
    end
end

y_coor_left_h_boundary_sorted = sort(y_coor(left_h_boundary_nodes));
y_coor_right_h_boundary_sorted = sort(y_coor(right_h_boundary_nodes));
y_coor_top_h_boundary_sorted = sort(y_coor(top_h_boundary_nodes));
y_coor_bottom_h_boundary_sorted = sort(y_coor(bottom_h_boundary_nodes));
y_coor_circle_h_boundary_sorted = sort(y_coor(circle_h_boundary_nodes));

n_eq = counter;

LM = ID(IEN);

h_integration = zeros(n_np, 1);

% left h boundary
for ee = 1 : length(y_coor_left_h_boundary_sorted)-1
    h_ele = zeros(2,1);
    x_ele = [y_coor_left_h_boundary_sorted(ee), y_coor_left_h_boundary_sorted(ee+1)];
    for qua = 1 : n_int_h
        dx_dxi_h = 0.0;
        x_l_h = 0.0;
        for aa = 1 : 2
            x_l_h =    x_l_h    + x_ele(aa) * PolyShape(1, aa, xi_h(qua), 0);
            dx_dxi_h = dx_dxi_h + x_ele(aa) * PolyShape(1, aa, xi_h(qua), 1);
        end
        dxi_dx_h = 1.0 / dx_dxi_h;
        for aa = 1 : 2
            h_ele(aa) = h_ele(aa) + weight_h(qua) * PolyShape(1, aa, xi_h(qua), 0) * h_left(x_l_h) * dx_dxi_h;
        end
    end
    for aa = 1 : 2
        for i = 1 : n_np
            if x_coor(i) == 1  && y_coor(i) == x_ele(aa)
                A = i;
                break;
            end
        end
        h_integration(A) = h_integration(A) + h_ele(aa);
    end
end

% right h boundary
for ee = 1 : length(y_coor_right_h_boundary_sorted)-1
    h_ele = zeros(2,1);
    x_ele = [y_coor_right_h_boundary_sorted(ee), y_coor_right_h_boundary_sorted(ee+1)];
    for qua = 1 : n_int_h
        dx_dxi_h = 0.0;
        x_l_h = 0.0;
        for aa = 1 : 2
            x_l_h =    x_l_h    + x_ele(aa) * PolyShape(1, aa, xi_h(qua), 0);
            dx_dxi_h = dx_dxi_h + x_ele(aa) * PolyShape(1, aa, xi_h(qua), 1);
        end
        dxi_dx_h = 1.0 / dx_dxi_h;
        for aa = 1 : 2
            h_ele(aa) = h_ele(aa) + weight_h(qua) * PolyShape(1, aa, xi_h(qua), 0) * h_right(x_l_h) * dx_dxi_h;
        end
    end
    for aa = 1 : 2
        for i = 1 : n_np
            if x_coor(i) == 1  && y_coor(i) == x_ele(aa)
                A = i;
                break;
            end
        end
        h_integration(A) = h_integration(A) + h_ele(aa);
    end
end

%top h boundary
for ee = 1 : length(y_coor_top_h_boundary_sorted)-1
    h_ele = zeros(2,1);
    x_ele = [y_coor_top_h_boundary_sorted(ee), y_coor_top_h_boundary_sorted(ee+1)];
    for qua = 1 : n_int_h
        dx_dxi_h = 0.0;
        x_l_h = 0.0;
        for aa = 1 : 2
            x_l_h =    x_l_h    + x_ele(aa) * PolyShape(1, aa, xi_h(qua), 0);
            dx_dxi_h = dx_dxi_h + x_ele(aa) * PolyShape(1, aa, xi_h(qua), 1);
        end
        dxi_dx_h = 1.0 / dx_dxi_h;
        for aa = 1 : 2
            h_ele(aa) = h_ele(aa) + weight_h(qua) * PolyShape(1, aa, xi_h(qua), 0) * h_top(x_l_h) * dx_dxi_h;
        end
    end
    for aa = 1 : 2
        for i = 1 : n_np
            if x_coor(i) == 1  && y_coor(i) == x_ele(aa)
                A = i;
                break;
            end
        end
        h_integration(A) = h_integration(A) + h_ele(aa);
    end
end

% bottom h boudary
for ee = 1 : length(y_coor_bottom_h_boundary_sorted)-1
    h_ele = zeros(2,1);
    x_ele = [y_coor_bottom_h_boundary_sorted(ee), y_coor_bottom_h_boundary_sorted(ee+1)];
    for qua = 1 : n_int_h
        dx_dxi_h = 0.0;
        x_l_h = 0.0;
        for aa = 1 : 2
            x_l_h =    x_l_h    + x_ele(aa) * PolyShape(1, aa, xi_h(qua), 0);
            dx_dxi_h = dx_dxi_h + x_ele(aa) * PolyShape(1, aa, xi_h(qua), 1);
        end
        dxi_dx_h = 1.0 / dx_dxi_h;
        for aa = 1 : 2
            h_ele(aa) = h_ele(aa) + weight_h(qua) * PolyShape(1, aa, xi_h(qua), 0) * h_bottom(x_l_h) * dx_dxi_h;
        end
    end
    for aa = 1 : 2
        for i = 1 : n_np
            if x_coor(i) == 1  && y_coor(i) == x_ele(aa)
                A = i;
                break;
            end
        end
        h_integration(A) = h_integration(A) + h_ele(aa);
    end
end

% circle h boudary
for ee = 1 : length(y_coor_circle_h_boundary_sorted)-1
    h_ele = zeros(2,1);
    x_ele = [y_coor_circle_h_boundary_sorted(ee), y_coor_circle_h_boundary_sorted(ee+1)];
    for qua = 1 : n_int_h
        dx_dxi_h = 0.0;
        x_l_h = 0.0;
        for aa = 1 : 2
            x_l_h =    x_l_h    + x_ele(aa) * PolyShape(1, aa, xi_h(qua), 0);
            dx_dxi_h = dx_dxi_h + x_ele(aa) * PolyShape(1, aa, xi_h(qua), 1);
        end
        dxi_dx_h = 1.0 / dx_dxi_h;
        for aa = 1 : 2
            h_ele(aa) = h_ele(aa) + weight_h(qua) * PolyShape(1, aa, xi_h(qua), 0) * h_circle(x_l_h) * dx_dxi_h;
        end
    end
    for aa = 1 : 2
        for i = 1 : n_np
            if x_coor(i) == 1  && y_coor(i) == x_ele(aa)
                A = i;
                break;
            end
        end
        h_integration(A) = h_integration(A) + h_ele(aa);
    end
end

% allocate the stiffness matrix and load vector
%K = spalloc(n_eq, n_eq, 9 * n_eq);
K = zeros(n_eq, n_eq);
F = zeros(n_eq, 1);

% loop over element to assembly the matrix and vector
for ee = 1 : n_el
    x_ele = x_coor( IEN(ee, 1:n_en) );    % 局部节点的全局横坐标 = x_ele(a)
    y_ele = y_coor( IEN(ee, 1:n_en) );    % 局部节点的全局纵坐标 = y_ele(a)

    k_ele = zeros(n_en, n_en); % element stiffness matrix
    f_ele = zeros(n_en, 1);    % element load vector
    g_ele = zeros(n_en, 1);    % element g boundary contition

    for ll = 1 : n_int
        x_l = 0.0; y_l = 0.0;
        dx_dxi = 0.0; dx_deta = 0.0;
        dy_dxi = 0.0; dy_deta = 0.0;
        for aa = 1 : n_en
            x_l = x_l + x_ele(aa) * Tria(aa, xi(ll), eta(ll));  %坐标变换Mapping
            y_l = y_l + y_ele(aa) * Tria(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Tria_grad(aa, xi(ll), eta(ll));
            dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
            dx_deta = dx_deta + x_ele(aa) * Na_eta;
            dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
            dy_deta = dy_deta + y_ele(aa) * Na_eta;
        end

        detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;

        for aa = 1 : n_en
            Na = Tria(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Tria_grad(aa, xi(ll), eta(ll));
            Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
            Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;

            f_ele(aa) = f_ele(aa) + weight(ll) * detJ * f(x_l, y_l) * Na;

            for bb = 1 : n_en
                Nb = Tria(bb, xi(ll), eta(ll));
                [Nb_xi, Nb_eta] = Tria_grad(bb, xi(ll), eta(ll));
                Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
                Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;

                k_ele(aa, bb) = k_ele(aa,bb) + weight(ll) * detJ * kappa * (Na_x * Nb_x + Na_y * Nb_y);
                if LM(ee, bb) == 0
                    g_ele(aa) = g_ele(aa) + weight(ll) * detJ * kappa * (Na_x * Nb_x + Na_y * Nb_y)*gg(IEN(ee,bb));
                end

            end % end of bb loop
        end % end of aa loop
    end % end of quadrature loop

    for aa = 1 : n_en
        PP = LM(ee, aa);
        if PP > 0
            F(PP) = F(PP) + f_ele(aa)- g_ele(aa);

            for bb = 1 : n_en
                QQ = LM(ee, bb);
                if QQ > 0
                    K(PP, QQ) = K(PP, QQ) + k_ele(aa, bb);

                else
                    % modify F with the boundary data
                    % here we do nothing because the boundary data g is zero or
                    % homogeneous
                end
            end
        end
    end
end

for AA = 1 : n_np
        PP = ID(AA);
        if PP > 0
            F(PP) = F(PP) + h_integration(AA);
        end
end

% solve the stiffness matrix
dn = K \ F;

% insert dn back into the vector for all nodes
disp = zeros(n_np, 1);

for ii = 1 : n_np
    index = ID(ii);
    if index > 0
        disp(ii) = dn(index);
    else
        % modify disp with the g data. Here it does nothing because g is zero
        disp(ii) = disp(ii) + gg(ii);
    end
end

% save the solution vector and number of elements to disp with name
% HEAT.mat

save("Solid-tria", "disp", "x_coor", "y_coor");

% error estimate
% e_0 = [];
% e_1 = [];
% for nel = 2:2:64
%     [e_0_temp, e_1_temp] = Erro_estimate_tria(nel);
%     e_0 = [e_0,e_0_temp];
%     e_1 = [e_1,e_1_temp];
% end
% nel = 2:2:64;
% nel = 1./nel;
% figure
% plot(log(nel),log(e_0),'LineWidth',1)
% xlabel("log(h_e_l)")
% ylabel("log(e_0\_error)")
% p=polyfit(log(nel),log(e_0),1);
% slope_e_0 = p(1)
% figure
% plot(log(nel),log(e_1),'LineWidth',1)
% xlabel("log(h_e_l)")
% ylabel("log(e_1\_error)")
% p=polyfit(log(nel),log(e_1),1);
% slope_e_1 = p(1)