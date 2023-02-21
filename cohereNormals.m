function [untreatedNormals,coheredNormals,residuals] = cohereNormals(points,r,c)
%cohereNormals Cohere orientations of normals using algorithm from http://hhoppe.com/recon.pdf
%   points is n x 3 array of x,y,z coords of points
%   y is num of nearest neighbors used to calc local normal at each point
%   c is num of those neighbors to connect the first graph

%parameters
%r = 10; %must be scaled already
min_num_neighbors = 3;


%aux
n = size(points,1);
normals = []; %indices match 'points' array, n x 3
residuals = [];
neighbors = []; % contains indices of each point's neighbor, indices match 'points', n x y

r_sq = r^2;

%get normals
for i = 1:n
    coord = points(i,:);
    
    %find cells within radius r of cell i
    neighbs = [];
    for j = setdiff([1:n],i)
        if sum((points(j,:)-coord).^2) < r_sq
            neighbs = [neighbs;j];
        end
    end
    %if not enough points in bounding sphere, take min_num_neighbors closest nns
    if length(neighbs) < min_num_neighbors
        [~,inds] = pdist2(points,coord,'euclidean','Smallest',3+1);
        neighbs=inds(2:end);
    end
    %log min_num_neighbors closest neighbs in 'neighbors'
    neighbors = [neighbors;neighbs(1:3)'];
    
    
    %local normal 
    [U,S,V]=svd(points(neighbs,:)','econ'); %note the transpose
    normal=U(:,1); % 3 x 1, unit normal %%****** changed to 1 because it works -- mar 22
    normals = [normals;normal'];
    
    residual = 0;
    for j = 1:length(neighbs)
        residual = residual + abs(dot(points(neighbs(j),:) - coord,normal));
    end
    residuals = [residuals;residual];
end
untreatedNormals = normals;



%troubleshooting
%frac_before_coherence = fracPointingOut(points,normals)



%get index of point with highest z value 
zval = 0;
zind = 0;
for i = 1:n
    if points(i,3) > zval
        zval = points(i,3);
        zind = i;
    end
end

%force orient the highest point's normal up
if normals(zind,3) < 0
    normals(zind,:) = -normals(zind,:);
end

%adjacency matrix
    %each node is connected to its c closest nns and all nodes for
    %which it is one of their c closest nns
    %weights are 1-dot(points' normals)
Adj = zeros(n);
for pt = 1:n
    for i = 1:c
        dex = neighbors(pt,i); %index of point pt's i'th closest nearest neighbor
        Adj(i,dex) = 1 - dot(normals(pt),normals(dex));
        Adj(dex,i) = 1 - dot(normals(pt),normals(dex));
    end
end

%generate path of least curvature through embryo
G = graph(Adj);
T = minspantree(G);
minAdj = adjacency(T);
minAdj_full = full(minAdj);
[disc,pred,closed] = graphtraverse(minAdj,zind);

%propogate normal directions through path
for i = 2:length(disc)
    %locate parent node
    backstep = 1;
    while minAdj_full(disc(i),disc(i-backstep)) == 0
        backstep = backstep + 1;
    end
    
    %make coherent
    if dot(normals(i-backstep,:),normals(i,:)) < 0
        normals(i,:) = -normals(i,:);
    end
end




%{
%for debugging
figure(1)
plot(G)
figure(2)
plot(T)
%}

%troubleshooting
%frac_after_coherence = fracPointingOut(points,normals)
coheredNormals = normals;

end

