clear all; clc;

kappa = 1.0; % conductivity
E = 1e9; % Young's modulus
upsilon = 0.3; % Poisson's ration υ
lambda = upsilon * E / ((1 + upsilon) * (1 - 2 * upsilon)); % λ
miu = E / (2 * (1 + upsilon)); % μ

n_sd = 2; % 固体力学问题的维度
n_dof = 2; % 节点自由度
D = [lambda+2*miu   lambda         0;
    lambda         lambda+2*miu   0;
    0              0              miu];
e = [1 0;
    0 1];  % e_1 = e(:,1)   e_2 = e(:, 2)




% exact solution
exact = @(x,y,i) (i == 1) * (x*(x-1)*y*(y-1)) + (i == 2) * (x*(x-1)*y*(y-1));
nu = 0.3;
f = @(x,y,i) (i == 1) * -( (E / (1 - nu^2)) * (2*y*(y - 1) + nu*( (x - 1)*(2*y - 1) + x*(2*y - 1) ) ) + ...
                           (E / (2*(1 + nu))) * ( x*(x - 1)*2 + (2*x - 1)*(x - 1) + 2*y*(2*x - 1) ) ) + ...
             (i == 2) * -( (E / (2*(1 + nu))) * ( (x - 1)*(2*y - 1) + y*(2*x - 1) + 2*y*(2*x - 1) )   + ...
                           (E / (1 - nu^2)) * (2*x*(x - 1) + nu*( (y - 1)*(2*x - 1) + y*(2*x - 1) ) ) );
Tx = 1e4;
R = 0.5;
xigema_r_r         = @(r,theta) Tx/2*(1-R^2/r^2)+Tx/2*(1-4*R^2/r^2+3*R^4/r^4)*cos(2*theta);
xigema_theta_theta = @(r,theta) Tx/2*(1+R^2/r^2)-Tx/2*(1          +3*R^4/r^4)*cos(2*theta);
xigema_r_theta     = @(r,theta) -Tx/2*(1+2*R^2/r^2-3*R^4/r^4)*sin(2*theta);

% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);


% quadrature rule for h boundary
n_int_h = 10;
[xi_h, weight_h] = Gauss(n_int_h, -1, 1);

% mesh数据导入
fileID = fopen("..\..\gmsh-files\manufactured-solution-quad.msh", 'r');
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


fileID = fopen("..\..\gmsh-files\manufactured-solution-quad.msh", 'r');

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
        if numRows == 9
            % if rowDataPadded(5) == 1
            %     temp_ien = rowDataPadded(end-3 : end)';
            %     sj = temp_ien(4);
            %     temp_ien(4) = temp_ien(2);
            %     temp_ien(2) = sj;
            %     IEN = [IEN; temp_ien];
            %     n_el = n_el + 1;
            % else
                IEN = [IEN; rowDataPadded(end-3 : end)'];
                n_el = n_el + 1;
            % end
        end
    end
    tline = fgetl(fileID);
end


fclose(fileID);
n_en   = 4;               % number of nodes in an element
n_np   = length(x_coor); % total number of nodal points


gg = zeros(n_np,1); % g boundary condition
left_boundary_nodes = [];
right_boundary_nodes = [];
bottom_boundary_nodes = [];
top_boundary_nodes = [];
circle_boundary_nodes = [];


% ID array
ID = zeros(n_sd, n_np);
counter = 0;
gg = zeros(n_sd, n_np);
hh = zeros(n_sd, n_np);

for nn = 1 : n_np
    if nn == 1
        ID(1, nn) = 0;
        gg(1, nn) = 0;
        ID(2, nn) = 0;
        gg(2, nn) = 0;

    elseif nn == 2
        ID(1, nn) = 0;
        gg(1, nn) = 0;
        ID(2, nn) = 0;
        gg(2, nn) = 0;

    elseif nn == 3
        ID(1, nn) = 0;
        gg(1, nn) = 0;
        ID(2, nn) = 0;
        gg(2, nn) = 0;

    elseif nn == 4
        ID(1, nn) = 0;
        gg(1, nn) = 0;
        ID(2, nn) = 0;
        gg(2, nn) = 0;

    elseif x_coor(nn) == 0  % left boundary, symmetric
        ID(1, nn) = 0;
        gg(1, nn) = 0;
        ID(2, nn) = 0;
        hh(2, nn) = 0;

    elseif y_coor(nn) == 0 % bottom boundary, symmetric
       ID(1, nn) = 0;
        gg(1, nn) = 0;
        ID(2, nn) = 0;
        hh(2, nn) = 0;

    elseif x_coor(nn) == 1  % right boundary, user define
        ID(1, nn) = 0;
        gg(1, nn) = 0;
        ID(2, nn) = 0;
        hh(2, nn) = 0;

    elseif y_coor(nn) == 1  % top boundary, user define
        ID(1, nn) = 0;
        gg(1, nn) = 0;
        ID(2, nn) = 0;
        hh(2, nn) = 0;

    else
        for ii = 1 : n_sd
            counter = counter +1;
            ID(ii, nn) = counter;
        end

    end
end

% right_boundary_nodes = [right_boundary_nodes; 1];
% top_boundary_nodes = [top_boundary_nodes; 2; 3];
% bottom_boundary_nodes = [bottom_boundary_nodes; 4];
% 
% h_right= @(y) [1e4 0];
% for ii = 1 : length(right_boundary_nodes)
%     hh(:, right_boundary_nodes(ii)) = hh(:, right_boundary_nodes(ii)) + h_right(y_coor(right_boundary_nodes(ii)))';
% end
% 
% 
% y_coor_right_h_boundary_sorted = sort(y_coor(right_boundary_nodes));
% 
% 
% 
n_eq = counter;
% 
% 
% h_integration = zeros(n_sd, n_np);
% 
% 
% % right h boundary
% for ee = 1 : length(y_coor_right_h_boundary_sorted)-1
%     for ii = 1 : 2
%         h_ele = zeros(2,1);
%         x_ele = [y_coor_right_h_boundary_sorted(ee), y_coor_right_h_boundary_sorted(ee+1)];
%         for qua = 1 : n_int_h
%             dx_dxi_h = 0.0;
%             x_l_h = 0.0;
%             for aa = 1 : 2
%                 x_l_h =    x_l_h    + x_ele(aa) * PolyShape(1, aa, xi_h(qua), 0);
%                 dx_dxi_h = dx_dxi_h + x_ele(aa) * PolyShape(1, aa, xi_h(qua), 1);
%             end
%             dxi_dx_h = 1.0 / dx_dxi_h;
%             for aa = 1 : 2
%                 h_temp = h_right(x_l_h);
%                 h_ele(aa) = h_ele(aa) + weight_h(qua) * PolyShape(1, aa, xi_h(qua), 0) * h_temp(ii) * dx_dxi_h;
%             end
%         end
%         for aa = 1 : 2
%             for i = 1 : n_np
%                 if x_coor(i) == 1  && y_coor(i) == x_ele(aa)
%                     A = i;
% 
%                     break;
%                 end
%             end
%             h_integration(ii, A) = h_integration(ii, A) + h_ele(aa);
%         end
%     end
% end


% allocate the stiffness matrix and load vector
%K = spalloc(n_eq, n_eq, 9 * n_eq);
K = zeros(n_eq, n_eq);
F = zeros(n_eq, 1);

% loop over element to assembly the matrix and vector
for ee = 1 : n_el
    for ii = 1 : n_sd
        for jj = 1 : n_sd
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
                    x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));  %坐标变换Mapping
                    y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
                    [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
                    dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
                    dx_deta = dx_deta + x_ele(aa) * Na_eta;
                    dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
                    dy_deta = dy_deta + y_ele(aa) * Na_eta;
                end

                detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;

                for aa = 1 : n_en
                    Na = Quad(aa, xi(ll), eta(ll));
                    [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
                    Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
                    Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
                    B_a = [Na_x 0;
                        0    Na_y;
                        Na_y Na_x];


                    f_ele(aa) = f_ele(aa) + weight(ll) * detJ * f(x_l, y_l, ii) * Na;
                    for bb = 1 : n_en
                        Nb = Quad(bb, xi(ll), eta(ll));
                        [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
                        Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
                        Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;
                        B_b = [Nb_x 0;
                            0    Nb_y;
                            Nb_y Nb_x];

                        k_ele(aa, bb) = k_ele(aa,bb) + weight(ll) * detJ * e(:, ii)' * B_a' * D * B_b * e(:, jj);
                        if ID(jj ,IEN(ee, bb)) == 0
                            g_ele(aa) = g_ele(aa) + weight(ll) * detJ * e(:, ii)' * B_a' * D * B_b * e(:, jj) * gg(jj, IEN(ee,bb));
                        end

                    end % end of bb loop
                end % end of aa loop
            end % end of quadrature loop

            for aa = 1 : n_en
                PP = ID(ii, IEN(ee, aa));
                if PP > 0
                    F(PP) = F(PP) + f_ele(aa) - g_ele(aa);

                    for bb = 1 : n_en
                        QQ = ID(jj, IEN(ee, bb));
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
    end
end

% for AA = 1 : n_np
%     for ii = 1 : n_sd
%         PP = ID(ii, AA);
%         if PP > 0
%             F(PP) = F(PP) + h_integration(ii, AA);
%         end
%     end
% end

% solve the stiffness matrix
dn = K \ F;

% insert dn back into the vector for all nodes
disp = zeros(n_sd, n_np);

for AA = 1 : n_np
    for ii = 1 : n_sd
        index = ID(ii, AA);
        if index > 0
            disp(ii, AA) = dn(index);
        else
            % modify disp with the g data. Here it does nothing because g is zero
            disp(ii, AA) = disp(ii, AA) + gg(ii, AA);
        end
    end
end

% save the solution vector and number of elements to disp with name
% HEAT.mat
save("Solid-quad", "disp", "x_coor", "y_coor");

% EOF

%数据可视化
% n_sam = 10; %每个单元绘图样本点的个数
% xi_sam = -1 : (2/n_sam) : 1;
% yi_sam = -1 : (2/n_sam) : 1;
%
% x_sam = zeros(n_el * n_sam + 1, 1);
% y_sam = zeros(n_el * n_sam + 1, 1);
% u_sam = zeros(n_el * n_sam + 1, n_el * n_sam + 1); % store the exact solution value at sampling points
% uh_sam = zeros(n_el * n_sam + 1, n_el * n_sam + 1); % store the numerical solution value at sampling pts
%
% for ex = 1 : n_el_x
%     for ey = 1 : n_el_y
%         ee = (ey-1) * n_el_x + ex;
%         x_ele = x_coor( IEN(ee, :) );
%         y_ele = y_coor( IEN(ee, :) );
%         u_ele = disp( IEN(ee, :) );
%
%         if ex == n_el_x % 最后一个单元多绘制一个点
%             n_sam_end_x = n_sam+1;
%         else
%             n_sam_end_x = n_sam;
%         end
%
%         if ey == n_el_y % 最后一个单元多绘制一个点
%             n_sam_end_y = n_sam+1;
%         else
%             n_sam_end_y = n_sam;
%         end
%
%         for ll = 1 : n_sam_end_x
%             for kk = 1 : n_sam_end_y
%                 x_l = 0.0;
%                 y_l = 0.0;
%                 u_l = 0.0;
%                 for aa = 1 : n_en
%                     x_l = x_l + x_ele(aa) * Quad( aa, xi_sam(ll),yi_sam(kk)); % 局部向全局的坐标变换
%                     y_l = y_l + y_ele(aa) * Quad( aa, xi_sam(ll),yi_sam(kk)); % 局部向全局的坐标变换
%                     u_l = u_l + u_ele(aa) * Quad( aa, xi_sam(ll),yi_sam(kk)); % u(x)解的表达式
%                 end
%
%                 x_sam( (ex-1)*n_sam + ll ) = x_l;
%                 y_sam( (ey-1)*n_sam + kk ) = y_l;
%                 uh_sam( (ex-1)*n_sam + ll, (ey-1)*n_sam + kk ) = u_l;
%                 u_sam( (ex-1)*n_sam + ll, (ey-1)*n_sam + kk ) = exact(x_l,y_l);
%             end
%         end
%     end
% end
% figure
% surf(x_sam, y_sam, u_sam);
% title("exact solution")
% xlabel("x")
% ylabel("y")
% zlabel("u")
% shading interp
% figure
% surf(x_sam, y_sam, uh_sam);
% title("Quad numercial solution")
% xlabel("x")
% ylabel("y")
% zlabel("u^h")
% shading interp

% error estimate
% e_0 = [];
% e_1 = [];
% for nel = 2:2:64
%     [e_0_temp, e_1_temp] = Erro_estimate_quia(nel);
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


disp_exact = zeros(2, n_np);
for AA = 1 : n_np
    for ii = 1 : 2
        disp_exact(ii, AA) = exact(x_coor(AA), y_coor(AA), ii);
    end
end