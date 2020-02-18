function [ H ] = est_homography(video_pts, logo_pts)
% est_homography estimates the homography to transform each of the
% video_pts into the logo_pts
% Inputs:
%     video_pts: a 4x2 matrix of corner points in the video
%     logo_pts: a 4x2 matrix of logo points that correspond to video_pts
% Outputs:
%     H: a 3x3 homography matrix such that logo_pts ~ H*video_pts
% Written for the University of Pennsylvania's Robotics:Perception course

% YOUR CODE HERE
H = [];
A = [];
for i = 1:length(video_pts)
    ax = [-video_pts(i,1) -video_pts(i,2) -1 0 0 0 video_pts(i,1)*logo_pts(i,1) ...
        video_pts(i,2)*logo_pts(i,1) logo_pts(i,1)];
    ay = [0 0 0 -video_pts(i,1) -video_pts(i,2) -1 video_pts(i,1)*logo_pts(i,2) ...
        video_pts(i,2)*logo_pts(i,2) logo_pts(i,2)];
    A = [A;ax;ay];
end
[U, S, V] = svd(A);
h = V(:,size(V,2));
H = [h(1:3,1)';h(4:6,1)';h(7:9,1)'];

end

