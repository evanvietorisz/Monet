function [ featurevector ] = computeFeaturesTimepoint(points,k)
%compute feature vector on each point in individual point cloud;
%note: points is n x 3 where n is determined by the state of the embryo
%passed in
featurevector=[];
for i=1:size(points,1)
    
    % distance to k nns
    [distfeat,inds]=pdist2(points,points(i,:),'euclidean','Smallest',k+1);
    distfeat=distfeat(2:end,:); % k x 1
    inds=inds(2:end); % k x 1
    
    %local normal 
    [coeff,score,roots]=svd(points(inds,:)','econ'); %note the transpose
    normal=coeff(:,3); % 3 x 1, unit normal
    
    %vectors pointing from anchor to k-1 nns
    directions=points(inds,:)-repmat(points(i,:),k,1); % k x 3
    
    %normalization
    %directions=directions./repmat(((directions(:,1).^2+directions(:,2).^2+directions(:,3).^2).^.5),1,3);
    
    %projection of directions onto normal
    normscalproj = dot(directions,repmat(normal',k,1),2); % k x 1
    
    normvecproj = normscalproj.*repmat(normal',k,1); % k x 3

    normprojfeat = [];
    for i = 1:k
        normprojfeat = [normprojfeat;normvecproj(i,:)']; % makes normprojfeat k*3 x 1
    end

    
    %projection of directions onto plane orthogonal to normal
    perpvecproj = directions-normvecproj; % k x 3
    perpscalproj = dot(perpvecproj,perpvecproj,2).^(1/2); % k x 1
    
    perpprojfeat = [];
    for i = 1:k
        perpprojfeat = [perpprojfeat;perpvecproj(i,:)']; %makes perpprojfeat k*3 x 1
    end
    
    
    %cos distance between anchor and k nns
    distcos=pdist2(directions,normal','cosine'); % k x 1

    %cos distance between anchor k nns
    %distcos=pdist2(points(inds,:),points(i,:),'cosine')';
        % A doesn't think it's a useful feature
        
        
    %{    
    %thrift-like feature
    k_small = k/2;
    
    [small_distfeat,small_inds]=pdist2(points,points(i,:),'euclidean','Smallest',k_small+1);
    small_distfeat=small_distfeat(2:end,:); % k_small x 1
    small_inds=small_inds(2:end); % k_small x 1

    
    %local(er) normal 
    [U,S,V]=svd(points(small_inds,:)','econ'); %note the transpose
    small_normal=U(:,3); % 3 x 1, unit normal
    
    thrift = dot(small_normal,normal); %both already normalized
    
    [thriftfeat, edges] = histcounts(thrift,[-1 -.5 0 .5 1]);
    %}
    

    % note: distfeat, inds,and distcos are k x 1.
        % distfeat is the euclidean distances to k
        % nearest neighbors.
        % distcos is cosine distances to k nns where
        % the cosine distance is 1 - cos(theta) where 
        % theta is the angle between the vector going from
        % the nearest neighbor to the point i in question
        % and the local normal calculated from the svd. 

    featurevector=[featurevector;distfeat',distcos'];
    %featurevector=[featurevector;distfeat',distcos',normprojfeat',perpprojfeat']; %take transpose
   
    %note: In featurevector, diff rows are diff points and columns are
        % 'blocks' of information k columns long that show (in this case)
        % distances to k nearest neighbors for a certain point.
        % Interpreting the feature data meaningfully means keeping track of
        % how many columns are used for each feature

end

end

