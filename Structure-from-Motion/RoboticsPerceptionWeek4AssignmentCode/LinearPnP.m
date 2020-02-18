function [C, R] = LinearPnP(X, x, K)
%% LinearPnP
% Getting pose from 2D-3D correspondences
% Inputs:
%     X - size (N x 3) matrix of 3D points
%     x - size (N x 2) matrix of 2D points whose rows correspond with X
%     K - size (3 x 3) camera calibration (intrinsics) matrix
% Outputs:
%     C - size (3 x 1) pose transation
%     R - size (3 x 3) pose rotation
%
% IMPORTANT NOTE: While theoretically you can use the x directly when solving
% for the P = [R t] matrix then use the K matrix to correct the error, this is
% more numeically unstable, and thus it is better to calibrate the x values
% before the computation of P then extract R and t directly

N = size(X,1);
x = [x ones(N,1)];
x = x';
x = K\x;
x = x';
A = [];


for i = 1:N
    % concatenate all equations together and solve with SVD
    X1 = X(i,1);
    X2 = X(i,2);
    X3 = X(i,3);
    x1 = x(i,1);
    x2 = x(i,2);
    ax = [-X1,-X2,-X3,0,0,0,X1*x1,X2*x1,X3*x1,-1,0,x1];
    ay = [0,0,0,-X1,-X2,-X3,X1*x2,X2*x2,X3*x2,0,-1,x2];
    A = [A;ax;ay];
end

[~,~,V] = svd(A);
params = V(:,size(V,2));

% recover R and t
R = reshape(params(1:end-3),[3,3])';
t = params(end-2:end);

[U,D,V] = svd(R);
if abs(det(U*V')-1) < 1e-3
    R = U*V';
    t = t/D(1,1);
elseif abs(det(U*V')+1) < 1e-3
    R = -U*V';
    t = -t/D(1,1);
end

C = -R'*t;

