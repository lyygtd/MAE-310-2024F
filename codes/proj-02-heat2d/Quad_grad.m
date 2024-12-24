function [val_xi, val_eta] = Quad_grad(aa, xi, eta)%型函数的导数
%返回值分别为 对ξ求导 和 对η求导 的结果
if aa == 1
    val_xi  = -0.25 * (1-eta);
    val_eta = -0.25 * (1-xi);
elseif aa == 2
    val_xi  =  0.25 * (1-eta);
    val_eta = -0.25 * (1+xi);
elseif aa == 3
    val_xi  = 0.25 * (1+eta);
    val_eta = 0.25 * (1+xi);
elseif aa == 4
    val_xi  = -0.25 * (1+eta);
    val_eta =  0.25 * (1-xi);
else
    error('Error: value of a should be 1,2,3, or 4.');
end

% EOF