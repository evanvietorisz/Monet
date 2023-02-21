
%load points
time = 305;
points = embdat_stabilized(time).finalpoints;
num_points = size(points,1);
k = 20; %number of nearest neighbors used to calculate normals
c = 5; %number of neighbors aligned to a given point in question
num_passes = 3; %number of passes used to propogate alignment

figure
scatter3(points(:,1),points(:,2),points(:,3))

%get normals
normals = [];
for i = 1:num_points

    %nearest neighbors
    % distance to k-1 nns
    [distfeat,inds]=pdist2(points,points(i,:),'euclidean','Smallest',k);
    distfeat=distfeat(2:end,:); % k-1 x 1
    inds=inds(2:end); % k-1 x 1
    
    %local normal 
    [U,S,V]=svd(points(inds,:)','econ'); %note the transpose
    normal=U(:,3); % 3 x 1, unit normal
    normals = [normals;normal'];
end

%figure out how normals oriented
points_plus_norms = points + normals;
orientations = zeros(num_points,1);
for i = 1:num_points
    %if normal points inward, set value to 1
    if norm(points(i)) >= norm(points_plus_norms(i))
        orientations(i) = 1;
    end
end

%fraction of normals oriented out
frac = 0;
for i = 1:num_points
    frac = frac + orientations(i);
end
frac_out = frac/num_points

%naive way of propagating outward orientation

%get index of point with highest z value 
zval = 0;
zind = 0;
for i = 1:num_points
    if points(i,3) > zval
        zval = points(i,3);
        zind = i;
    end
end

%force orient the highest point's normal up
if norm(points(zind,:)) > norm(points(zind,:)+normals(zind,:))
    normals(zind,:) = -normals(zind,:);
end

points_to_hit = [zind];

for p = 1:num_passes
    next_batch = [];
    for i = 1:length(points_to_hit)
        % get indices of p closest neighbors
        [distances,indices]=pdist2(points,points(points_to_hit(i),:),'euclidean','Smallest',c+1);
        indices=indices(2:end); % p x 1
        
        %for each close neighbor, align normals
        for j = 1:length(indices)
            if dot(normals(points_to_hit(i),:),normals(indices(j),:)) < 0
                normals(indices(j)) = -normals(indices(j));
            end
        end
        
        %keep track of the those neighbors' indices to to hit *their*
        %neighbors on the next pass
        next_batch = [next_batch, indices'];
    end
    points_to_hit = [points_to_hit, next_batch];
end

%figure out how normals oriented AGAIN
points_plus_norms = points + normals;
orientations = zeros(num_points,1);
for i = 1:num_points
    %if normal points inward, set value to 1
    if norm(points(i)) >= norm(points_plus_norms(i))
        orientations(i) = 1;
    end
end

%fraction of normals oriented out
frac = 0;
for i = 1:num_points
    frac = frac + orientations(i);
end
frac_out_after_alignment = frac/num_points


%method doesn't work-->use graph propogation technique
