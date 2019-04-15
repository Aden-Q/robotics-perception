function X = Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1, x2, x3, X0)
%% Nonlinear_Triangulation
% Refining the poses of the cameras to get a better estimate of the points
% 3D position
% Inputs: 
%     K - size (3 x 3) camera calibration (intrinsics) matrix
%     x
% Outputs: 
%     X - size (N x 3) matrix of refined point 3D locations 

X = X0;
for i = 1:length(X)
    while true
        X_new = Single_Point_Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1(i,:), x2(i,:), x3(i,:), X(i,:)')';
        if norm(X_new-X(i,:)) < 1e-8
            break
        else
            % X(i,:)
            % X_new
            X(i,:) = X_new;
        end
    end
end

end

function X = Single_Point_Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1, x2, x3, X0)

b = [x1(1), x1(2), x2(1), x2(2), x3(1), x3(2)]';
uvw1 =K*R1*(X0-C1);
u1 = uvw1(1);
v1 = uvw1(2);
w1 = uvw1(3);
uvw2 =K*R2*(X0-C2);
u2 = uvw2(1);
v2 = uvw2(2);
w2 = uvw2(3);
uvw3 =K*R3*(X0-C3);
u3 = uvw3(1);
v3 = uvw3(2);
w3 = uvw3(3);
fx = [u1/w1 v1/w1 u2/w2 v2/w2 u3/w3 v3/w3]';
J1 = Jacobian_Triangulation(C1, R1, K, X0);
J2 = Jacobian_Triangulation(C2, R2, K, X0);
J3 = Jacobian_Triangulation(C3, R3, K, X0);
J = [J1 J2 J3]';
delta_X = (J'*J)\J'*(b-fx);
X = X0 + delta_X;

end

function J = Jacobian_Triangulation(C, R, K, X)
uvw = K*R*(X-C);
u = uvw(1);
v = uvw(2);
w = uvw(3);
d = K*R;
du = d(1,:);
dv = d(2,:);
dw = d(3,:);
df = [(w*du-u*dw)/w^2; (w*dv-v*dw)/w^2];
J = df';

end
