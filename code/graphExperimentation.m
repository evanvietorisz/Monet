time = 305;
points = embdat_stabilized(time).finalpoints; % n x 3

%rng(2)
%points = rand(10,3);

n = size(points,1);
y = 10; %number of nearest neighbors used to calculate normals
c = 10; %number of initial connections to per node. c <= y

%get normals
normals = []; %indices match 'points' array, n x 3
neighbors = []; % contains indices of each point's neighbor, indices match 'points', n x y
for pt = 1:n
    %nearest neighbors
    % distance to k-1 nns
    [distfeat,inds]=pdist2(points,points(pt,:),'euclidean','Smallest',y+1);
    distfeat=distfeat(2:end,:); % y x 1
    inds=inds(2:end); % y x 1
    neighbors = [neighbors;inds'];
    
    %local normal 
    [U,S,V]=svd(points(inds,:)','econ'); %note the transpose
    normal=U(:,3); % 3 x 1, unit normal
    normals = [normals;normal'];
end


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
    %each node is connected to its c closest nn's and all nodes for
    %which it is one of their c closest nn's
    %weights are 1-dot(points' normals)
Adj = zeros(n);
for pt = 1:n
    for i = 1:c
        dex = neighbors(pt,i); %index of point pt's i'th closest nearest neighbor
        Adj(pt,dex) = 1 - dot(normals(pt),normals(dex));
        Adj(dex,pt) = 1 - dot(normals(pt),normals(dex));
    end
end

G = graph(Adj);
figure(1)
plot(G)
T = minspantree(G);
figure(2)
plot(T)

minAdj = adjacency(T);
minAdj_full = full(minAdj);


[disc,pred,closed] = graphtraverse(minAdj,zind);

%propogate normal directions
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


