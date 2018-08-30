function out = imageStitching(file1, file2)

landscape_1 = file1;
landscape_2 = file2;

[l1_frames, l1_des] = vl_sift(landscape_1);
[l2_frames, l2_des] = vl_sift(landscape_2);

% get the differences between all descriptors from both images
distances = dist2(double(l1_des')/255, double(l2_des')/255);
[M, ~] = size(distances);

% new matrix to store the ratios of min / 2nd min
ratios = zeros(M, 3);
% sort the distances to later find the least
% store the y coords of these distances for plotting
[sorted_dist, y_dist] = sort(distances,2);

% loop to divide min and 2nd min ratios
% also store the coordinates of both images
for m = 1 : M
    min_1 = y_dist(m, 1);
    ratio = sorted_dist(m,1)/sorted_dist(m,2);
    ratios(m, 1) = ratio;
    ratios(m, 2) = m; % y value for im1
    ratios(m, 3) = min_1; % y value for im2
end

max_inliers = 0;
H_max_inliers = zeros(3,3);

% calculate minimum trials
S = ceil((log(1 - 0.90)) / (log(1 - (0.2)^size(ratios,1) - 0.001)));

for i = 1 : S
    
    % need to keep selecting 4 random matches
    matches = randperm(size(ratios,1), 4);
    match_1 = ratios(matches(1),:); %the whole row with ratio, y, y'
    match_2 = ratios(matches(2),:);
    match_3 = ratios(matches(3),:);
    match_4 = ratios(matches(4),:);

    % finding x,y coordinates of l1
    xy_r1 = l1_frames(:,match_1(2));
    x1 = xy_r1(1,1); y1 = xy_r1(2,1);
    xy_r2 = l1_frames(:,match_2(2));
    x2 = xy_r2(1,1); y2 = xy_r2(2,1);
    xy_r3 = l1_frames(:,match_3(2));
    x3 = xy_r3(1,1); y3 = xy_r3(2,1);
    xy_r4 = l1_frames(:,match_4(2));
    x4 = xy_r4(1,1); y4 = xy_r4(2,1);

    % finding x,y coordinates of l2
    xy_t1 = l2_frames(:, match_1(3));
    xp1 = xy_t1(1,1); yp1 = xy_t1(2,1);
    xy_t2 = l2_frames(:, match_2(3));
    xp2 = xy_t2(1,1); yp2 = xy_t2(2,1);
    xy_t3 = l2_frames(:, match_3(3));
    xp3 = xy_t3(1,1); yp3 = xy_t3(2,1);
    xy_t4 = l2_frames(:, match_4(3));
    xp4 = xy_t4(1,1); yp4 = xy_t4(2,1);

    % creating matrix A to compute homography H
    A = [x1, y1, 1, 0, 0, 0, -xp1*x1, -xp1*y1, -xp1; 
        0, 0, 0, x1, y1, 1, -yp1*x1, -yp1*y1, -yp1; 
        x2, y2, 1, 0, 0, 0, -xp2*x2, -xp2*y2, -xp2; 
        0, 0, 0, x2, y2, 1, -yp2*x2, -yp2*y2, -yp2;  
        x3, y3, 1, 0, 0, 0, -xp3*x3, -xp3*y3, -xp3; 
        0, 0, 0, x3, y3, 1, -yp3*x3, -yp3*y3, -yp3; 
        x4, y4, 1, 0, 0, 0, -xp4*x4, -xp4*y4, -xp4; 
        0, 0, 0, x4, y4, 1, -yp4*x4, -yp4*y4, -yp4;];
    
    % computing H
    [~, ~, V] = svd(A);
    h = double(vpa(V(:,9)));
    H = reshape(h,[3,3])';
    
    % checking for inliers by iterating through all matches
    inliers = 0;
    for j = 1 : size(ratios,1) % up till the number of rows, aka matches
        
        % get x,y of match in landscape_1
        match_l1 = ratios(j,:);
        xy_l1 = l1_frames(:,match_l1(2));
        x_l1 = xy_l1(1,1);
        y_l1 = xy_l1(2,1);
        
        % get x,y of match in landscape_2
        match_l2 = ratios(j,:);
        xy_l2 = l2_frames(:,match_l2(3));
        x_l2 = xy_l2(1,1);
        y_l2 = xy_l2(2,1);
        
        % map landscape_1 coordinate via H
        Hx_l1 = ((H(1,1)*x_l1 + H(1,2)*y_l1 + H(1,3))/ (H(3,1)*x_l1 + H(3,2)*y_l1 + H(3,3)));
        Hy_l1 = ((H(2,1)*x_l1 + H(2,2)*y_l1 + H(2,3))/ (H(3,1)*x_l1 + H(3,2)*y_l1 + H(3,3)));
        
        diff = pdist([x_l2, y_l2; Hx_l1, Hy_l1], 'euclidean');
        
        % increment inliers accordingly
        if (diff <= 3)
            inliers = inliers + 1;
        end
    end
    
    % compare with max_inliers
    if (inliers > max_inliers)
        max_inliers = inliers;
        H_max_inliers = H;
    end
end

% creating panaroma
% learned from http://home.deib.polimi.it/boracchi/teaching/IAS/Stitching/stitch.html
% maketform only used to convert my H matrix to type tform
T = maketform('projective',inv(H_max_inliers)');
[im2t,xdataim2t,ydataim2t] = imtransform(im,T);

% now xdataim2t and ydataim2t store the bounds of the transformed im2
xdataout=[min(1,xdataim2t(1)) max(size(im,2),xdataim2t(2))];
ydataout=[min(1,ydataim2t(1)) max(size(im,1),ydataim2t(2))];

% let's transform both images with the computed xdata and ydata
im2t=imtransform(im1,T, 'XData',xdataout,'YData',ydataout);
im1t=imtransform(im,maketform('affine',eye(3)), 'XData',xdataout,'YData',ydataout);
disp(size(im2t));
disp(size(im1t));

out=max(im1t,im2t);
imshow(out);

end