function [featurevector] = darbouxFeaturesTimepoint(points,k,kNorm,inputVector)
% Compute and return features as outlined in http://citeseerx.ist.psu.edu/viewdoc/
    %download;jsessionid=B0D59DA48B9CCC60EC9430B67D7C3B3B?doi=10.1.1.324.3396&rep=rep1&type=pdf

    %whereas above computes features for each *pair* of points, this does
    %it for every point
    
    %inputVector argument allows this to add to a featurevector computed in
    %another function
featurevector = [];
n = size(points,1);

%store indices of each point's nearest neighbors and each point's local normal
neighbor_inds = []; % k-1 x n, each column different point
[normals,treatedNormals] = cohereNormals(points,kNorm,5); % n x 3
normals = normals;
for i = 1:n
    % distance to k nns
    [dist_temp,inds]=pdist2(points,points(i,:),'euclidean','Smallest',k+1);
    dist_temp=dist_temp(2:end,:); % k x 1
    inds=inds(2:end); % k x 1
    
    neighbor_inds = [neighbor_inds,inds];
end

%Darboux
for i = 1:n
    %to avoid computing it twice
    neighbor_coords = points(neighbor_inds(:,i),:); % k x 3
    
    %vectors pointing from anchor to k nns
    directions=neighbor_coords-repmat(points(i,:),k,1); % k x 3
    
    %features
    feat_1 = []; % becomes k x 1
    feat_2 = []; % becomes k x 1
    feat_3 = []; % becomes k x 1
    feat_4 = []; % becomes k x 1
    for j = 1:k
        %darboux frame
        u = normals(i,:); % 1 x 3
        v = cross(directions(j,:),u);  % 1 x 3
        w = cross(u,v);  % 1 x 3
        
        val = dot(v,neighbor_coords(j,:));
        feat_1 = [feat_1;val];
        
        mag = norm(directions(j,:));
        feat_2 = [feat_2;mag];
        
        proj = dot(u,directions(j,:))/mag;
        feat_3 = [feat_3;proj];
        
        ang = atan(dot(w,normals(j,:)));
        feat_4 = [feat_4;ang];
    end
    
    featurevector = [featurevector;inputVector(i,:) feat_1',feat_2',feat_3',feat_4']; % n x (input + k*(# of features))

end
end

