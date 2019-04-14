function [proj_points, t, R] = ar_cube(H,render_points,K)
%% ar_cube
% Estimate your position and orientation with respect to a set of 4 points on the ground
% Inputs:
%    H - the computed homography from the corners in the image
%    render_points - size (N x 3) matrix of world points to project
%    K - size (3 x 3) calibration matrix for the camera
% Outputs: 
%    proj_points - size (N x 2) matrix of the projected points in pixel
%      coordinates
%    t - size (3 x 1) vector of the translation of the transformation
%    R - size (3 x 3) matrix of the rotation of the transformation
% Written by Stephen Phillips for the Coursera Robotics:Perception course

% YOUR CODE HERE: Extract the pose from the homography
if H(3,3) < 0
    H = -H;
end

h1 = H(:,1);
h2 = H(:,2);
h3 = H(:,3);
[U,~,V] = svd([h1,h2,cross(h1,h2)]);
temp = eye(3);
temp(3,3) = det(U*V');
R = U * temp * V';
t = h3 / norm(h1);

proj_points = zeros(size(render_points,1),2);

% YOUR CODE HERE: Project the points using the pose
for i=1:size(render_points,1)
    X = K*(R*render_points(i,:)'+t);
    x_im = X(1)/X(3);
    y_im = X(2)/X(3);
    proj_points(i,:) = [x_im,y_im];  
end

end
