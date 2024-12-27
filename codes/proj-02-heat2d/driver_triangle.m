clear all; clc;

kappa = 1.0; % conductivity

% exact solution
exact = @(x,y) x*(1-x)*y*(1-y);
exact_x = @(x,y) (1-2*x)*y*(1-y);
exact_y = @(x,y) x*(1-x)*(1-2*y);

f = @(x,y) 2.0*kappa*x*(1-x) + 2.0*kappa*y*(1-y); % source term

% quadrature rule
n_int = 3;
weight = [1/3, 1/3, 1/3]/2;
xi =     [2/3, 1/6, 1/6];
eta =    [1/6, 2/3, 1/6];


% mesh generation
n_en   = 3;               % number of nodes in an element
n_el_x = 3;               % number of elements in x-dir  直接将四边形网格按对角线一分为二的三角形网格,网格数量为原来的两倍
n_el_y = 3;               % number of elements in y-dir
n_el   = n_el_x * n_el_y * 2; % total number of elements

n_np_x = n_el_x + 1;      % number of nodal points in x-dir
n_np_y = n_el_y + 1;      % number of nodal points in y-dir
n_np   = n_np_x * n_np_y; % total number of nodal points

x_coor = zeros(n_np, 1);
y_coor = x_coor;

hx = 1.0 / n_el_x;        % mesh size in x-dir
hy = 1.0 / n_el_y;        % mesh size in y-dir

n_el_x = n_el_x * 2;
n_el_y = n_el_y * 2;

% generate the nodal coordinates
for ny = 1 : n_np_y
  for nx = 1 : n_np_x
    index = (ny-1)*n_np_x + nx; % nodal index
    x_coor(index) = (nx-1) * hx;
    y_coor(index) = (ny-1) * hy;
  end
end

% IEN array
IEN = zeros(n_el, n_en);
for ex = 1 : n_el_x
  for ey = 1 : n_el_y/2
    ee = (ey-1) * n_el_x + ex; % element index
    if mod(ey,2) == 1
        if mod(ex,4) == 1
            IEN(ee, 1) = (ey-1) * n_np_x + 1 + floor(ex / 4) * 2;
            IEN(ee, 2) =  ey    * n_np_x + 1 + floor(ex / 4) * 2;
            IEN(ee, 3) =  ey    * n_np_x + 2 + floor(ex / 4) * 2;
        end
        if mod(ex,4) == 2
            IEN(ee, 1) =  (ey-1) * n_np_x + 1 + floor(ex / 4) * 2;
            IEN(ee, 2) =  (ey-1) * n_np_x + 2 + floor(ex / 4) * 2;
            IEN(ee, 3) =   ey    * n_np_x + 2 + floor(ex / 4) * 2;
        end
        if mod(ex,4) == 3
            IEN(ee, 1) = (ey-1) * n_np_x + 2 + floor(ex / 4) * 2;
            IEN(ee, 2) = (ey-1) * n_np_x + 3 + floor(ex / 4) * 2;
            IEN(ee, 3) =  ey    * n_np_x + 2 + floor(ex / 4) * 2;
        end
        if mod(ex,4) == 0
            IEN(ee, 1) = (ey-1) * n_np_x + 3 + (ex / 4 - 1) * 2;
            IEN(ee, 2) =  ey    * n_np_x + 2 + (ex / 4 - 1) * 2;
            IEN(ee, 3) =  ey    * n_np_x + 3 + (ex / 4 - 1) * 2;
        end
    end
    if mod(ey,2) == 0
        if mod(ex,4) == 1
            IEN(ee, 1) = (ey-1) * n_np_x + 1 + floor(ex / 4) * 2;
            IEN(ee, 2) = (ey-1) * n_np_x + 2 + floor(ex / 4) * 2;
            IEN(ee, 3) =  ey    * n_np_x + 1 + floor(ex / 4) * 2;
        end
        if mod(ex,4) == 2
            IEN(ee, 1) =  (ey-1) * n_np_x + 2 + floor(ex / 4) * 2;
            IEN(ee, 2) =   ey    * n_np_x + 1 + floor(ex / 4) * 2;
            IEN(ee, 3) =   ey    * n_np_x + 2 + floor(ex / 4) * 2;
        end
        if mod(ex,4) == 3
            IEN(ee, 1) = (ey-1) * n_np_x + 2 + floor(ex / 4) * 2;
            IEN(ee, 2) =  ey    * n_np_x + 2 + floor(ex / 4) * 2;
            IEN(ee, 3) =  ey    * n_np_x + 3 + floor(ex / 4) * 2;
        end
        if mod(ex,4) == 0
            IEN(ee, 1) = (ey-1) * n_np_x + 2 + (ex / 4 - 1) * 2;
            IEN(ee, 2) = (ey-1) * n_np_x + 3 + (ex / 4 - 1) * 2;
            IEN(ee, 3) =  ey    * n_np_x + 3 + (ex / 4 - 1) * 2;
        end
    end
  end
end

% ID array
ID = zeros(n_np,1);
counter = 0;
for ny = 2 : n_np_y - 1
  for nx = 2 : n_np_x - 1
    index = (ny-1)*n_np_x + nx;
    counter = counter + 1;
    ID(index) = counter;  
  end
end

n_eq = counter;

LM = ID(IEN);

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
      end % end of bb loop
    end % end of aa loop
  end % end of quadrature loop
 
  for aa = 1 : n_en
    PP = LM(ee, aa);
    if PP > 0
      F(PP) = F(PP) + f_ele(aa);
      
      for bb = 1 : n_en
        QQ = LM(ee, bb);
        if QQ > 0
          K(PP, QQ) = K(PP, QQ) + k_ele(aa, bb);
          % fprintf("(PP,QQ)=(%d,%d)\n",PP,QQ)
          % fprintf("(ee,aa,bb)=(%d,%d,%d)\n",ee,aa,bb)
          % fprintf("k_ele(aa,bb)=%.1f\n",k_ele(aa,bb))
        else
          % modify F with the boundary data
          % here we do nothing because the boundary data g is zero or
          % homogeneous
        end
      end  
    end
  end
  % K
  % k_ele
  % fprintf('^^^^^^^^^^^^^^^^^^^^^^^\n')
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
  end
end

% save the solution vector and number of elements to disp with name
% HEAT.mat
save("HEAT", "disp", "n_el_x", "n_el_y");

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
%                     x_l = x_l + x_ele(aa) * Tria( aa, xi_sam(ll),yi_sam(kk)); % 局部向全局的坐标变换
%                     y_l = y_l + y_ele(aa) * Tria( aa, xi_sam(ll),yi_sam(kk)); % 局部向全局的坐标变换
%                     u_l = u_l + u_ele(aa) * Tria( aa, xi_sam(ll),yi_sam(kk)); % u(x)解的表达式
%                 end
% 
%                 x_sam( (ex-1)*n_sam + ll ) = x_l*hx;
%                 y_sam( (ey-1)*n_sam + kk ) = y_l*hy;
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
% figure
% surf(x_sam, y_sam, uh_sam);
% title("numercial solution")
% xlabel("x")
% ylabel("y")
% zlabel("u^h")
