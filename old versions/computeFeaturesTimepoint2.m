function [featurevector] = computeFeaturesTimepoint2(embryo,k,t)
%compute feature vector on each point in individual point cloud;

%points from embryo--non-temporal features work within one time step
points = embryo(t).finalpoints;

%num points
n = size(points,1);

%normals, calculated with k_norm nns
%***manually put 10 microns (conversion w xy scale); cohereNormals no
%longer uses kNorm nns
[normals,treatedNormals] = cohereNormals(points,10/.16,3);
norms = treatedNormals;

%compute centroid
cent = [0,0,0];
for i = n
    cent = cent + points(i,:);
end
cent = cent./n;

%time at which emb is 4 cells (developmental landmark)
four_cell_time = 1;
while size(embryo(four_cell_time).finalpoints,1) < 4
    four_cell_time = four_cell_time + 1;
end



%features
featurevector=[];
for i=1:n
    
    %if mod(t,10)==0 && i == 1
    if i==1
        ['time = ',num2str(t)]
    end
    
    %['time = ',num2str(t),', i = ',num2str(i)]
    %%% Local Features %%%
    %{
    % distance to k nns
    [distfeat,inds] = pdist2(points,points(i,:),'euclidean','Smallest',k+1);
    distfeat = distfeat(2:end,:); % k x 1
    inds = inds(2:end); % k x 1
    
    %vectors pointing from anchor to k nns
    directions = points(inds,:) - repmat(points(i,:),k,1); % k x 3
    
    %cos distance between anchor's normal and direction vectors to k nns
    distcos = pdist2(directions,norms(i,:),'cosine'); % k x 1
    
    %cos distance between anchor's normal and k nns' normals
    binorm = pdist2(norms(inds,:),norms(i,:),'cosine'); % k x 1
    
    %scalar projection of directions onto normal
    normscalproj = dot(directions,repmat(norms(i,:),k,1),2); % k x 1
    
    %vector projection of directions onto normal
    normvecproj = normscalproj.*repmat(norms(i,:),k,1); % k x 3
    %formal
    normprojfeat = [];
    for j = 1:k
        normprojfeat = [normprojfeat;normvecproj(j,:)']; % linearizes normvecproj k*3 x 1
    end
    
    %projection of directions onto plane orthogonal to normal
    perpvecproj = directions-normvecproj; % k x 3
    perpscalproj = dot(perpvecproj,perpvecproj,2).^(1/2); % k x 1
    
    perpprojfeat = [];
    for j = 1:k
        perpprojfeat = [perpprojfeat;perpvecproj(j,:)']; % linearizes perpvecproj k*3 x 1
    end
    
    %%% Global Features %%%
    
    %absolute position
    abs_position = points(i,:); % 1 x 3
    
    %relative position (absolute position - centroid)
    rel_position = points(i,:) - cent; % 1 x 3
    
    %%% Temporal Features %%%
    
    %absolute time
    time = t;
    %}
    
    %{
    %number of cells in emb
    num_cells = size(embryo(t).finalpoints,1);
    %}
    
    %{
    %migration path
    %   vector connecting ancestor x mins ago to cell in question
    %   actual sequence of paths as conceatenated horizontal vector
    num_mins_back = 30;
    migrationpath = [];
    
    steps_back = round(num_mins_back/1.26); %4cell->1.5 fold = 75mins->460mins; 385mins/305entries in embdat = 1.26 mins/entry
    
    current_time = t;
    current_cell_ind = i;
    
    for m = steps_back:-1:1
        for j = 1:size(embryo(current_time-1).suc,1)
            if embryo(current_time-1).suc(j,1) == current_cell_ind || embryo(current_time-1).suc(j,2) == current_cell_ind
                
                displacement = embryo(current_time).finalpoints(current_cell_ind,:) - embryo(current_time-1).finalpoints(j,:);
                migrationpath = [migrationpath,displacement];
                
                current_time = current_time - 1;
                current_cell_ind = j;
                break
            end
        end
    end
    
    displacementvec = embryo(t).finalpoints(i,:) - embryo(current_time).finalpoints(current_cell_ind,:);
    
    %directness ratio of piecewise migration path length to displacement
    pathlength = 0;
    for m = 1:3:length(migrationpath)
        pathlength = pathlength + norm(migrationpath(m:m+2));
    end
    directness = pathlength/norm(displacementvec);
   
   
    
    
    %{
    %developmental stage
    %0-gastrulation, gastrulation, gastrulation-elongation, elongation
    stage_edges = [0,26,330,400,460];
    %%%***need to figure out how to incorporate this binning with distance
    %%%features
    stage = histcounts(num_cells,stage_edges);
    %}
    
    %{
    %motion vector

    %   if a cell splits, motion vector points to avg position of children
    %   if cell dies, motion vector = [0,0,0]
    %   vector at a time t where each cell will go in the next timestep
    %   => not defined for last time step in dataset
    
  
    %grab last nonempty time
    nonemptime = length(embryo);
    while isempty(embryo(nonemptime).finalpoints)
        nonemptime = nonemptime - 1;
    end
    %if not penultimate meaningful time, compute
    if t < nonemptime
        if embryo(t).suc(i,1) == -1
            motion = [0,0,0];
        elseif embryo(t).suc(i,2) == -1
            nextpoint = embryo(t+1).finalpoints(embryo(t).suc(i,1),:);
            motion = nextpoint - points(i,:);
        else
            pt1 = embryo(t+1).finalpoints(embryo(t).suc(i,1),:);
            pt2 = embryo(t+1).finalpoints(embryo(t).suc(i,2),:);
            nextpoint = 1/2*(pt1 + pt2);
            motion = nextpoint - points(i,:);
        end
    %if final meaningful time, [0,0,0]
    else
        motion = [0,0,0];
    end
    %}
    %}
    
    
    
    
    
    
    %%
    %Cell Cycle
    
    %cell cycle metrics of current cell
    [beginning_of_current_cell_cycle,end_of_cell_current_cycle,duration_of_current_cell_cycle,time_until_next_division,time_since_last_division,fraction_of_current_cycle_elapsed,time_of_next_division] = cellCycleMetrics(t,i,embryo);
        %leave out time of next division; redundant
    

    
    
    
    
    %%

    %what is this for?
    %time of next division wrt 4 cell stage
    %***adjust for temp, UV light, same as developmental time
    scaled_time_of_next_division = time_of_next_division - four_cell_time;

    
    
    
    
    
    
    
    
    
    
    
    %%
    
    %Comparison with Cousins
    
    %hyperparameter
    %number of generations back to go (in this case both to find ancestors and cousins)
    g = 1;
    
    %traverse back g generations (cell divisions) and log the cell cycle metrics of
    %the ancestors along the way
    ancestor_metrics = [];
    
    time_pt = t;
    ind_pt = i;
    
    counter = g;
    skip_toggle = 0;
    while counter > 0
        %if you end up at beginning of the dataset, means that there are no
        %g'th cousins
        if time_pt == 1
            skip_toggle = 1;
            break
        end
        
        %walk back
        ind_pt = stepback(time_pt,ind_pt,embryo);
        time_pt = time_pt-1;
        %if crossed a cell division event, log information of progenitor,
        %increment counter
        if embryo(time_pt).suc(ind_pt,2) ~= -1
            [out1,out2,out3,out4,out5,out6,out7] = cellCycleMetrics(time_pt,ind_pt,embryo);
            ancestor_metrics = [ancestor_metrics;out1,out2,out3,out4,out5,out6,out7];
            counter = counter-1;
        end
        %leaves you off at time/index of progenitor cell g generations ago
        %at the last time it existed before dividing
        
    end
    
    
    if skip_toggle ~= 1 %if g'th cousins are defined, track forward
    
        %['time = ',num2str(t),', i = ',num2str(i)]
        %['time_pt = ',num2str(time_pt),' ind_pt = ',num2str(ind_pt)]
        
        %traverse forward back to the starting time having collected the
        %indices of the nth cousins
        cells_to_look_at = embryo(time_pt).suc(ind_pt,:); %guaranteed not to have -1 because of where time_pt,ind_pt point to (see above)
        time_pt = time_pt+1;
        
        
        [cells_to_look_at,time_pt] = forwardTraverse(time_pt,cells_to_look_at,'timesteps',t,embryo);
        
        %disp(find(cells_to_look_at==i))
        
        
        %remove index of cell i from cells_to_look_at
        cousin_indices = cells_to_look_at;
        cousin_indices(find(cousin_indices==i)) = [];
        
        
        %time_pt = time_pt+1; %leaves you off at time_pt = t****FIGURE OUT index off by one situation. not off by one here, but just double check
        %dont think this is necessary

        ind_pt = i; %time_pt,ind_pt now point back to t,i

        %compare cell i's cell cycle length to ancestors' and cousins'
        
    else %else, come back to t,i with no cousins
        time_pt = t;
        ind_pt = i;
        cousin_indices = [];     
    end
    
    
    
    %Comparison with ancestors
    [out1,out2,out3,out4,out5,out6,out7] = cellCycleMetrics(t,i,embryo);
    ancestor_metrics = [ancestor_metrics;out1,out2,out3,out4,out5,out6,out7]; %cell i is last row
    normed_ancestor_metrics = zscore(ancestor_metrics);
    similarity_to_ancestors = normed_ancestor_metrics(end,3); %hard-coded to grab duration of cell cycle column
    
    %cousins
    cousin_indices = [cousin_indices,i]; %i in last row
    %remove any dead cells that may have entered the cousins
    cousin_indices(find(cousin_indices==-1)) = [];
    
    cousin_metrics = [];
    for count = 1:length(cousin_indices)
        
        [out1,out2,out3,out4,out5,out6,out7] = cellCycleMetrics(t,cousin_indices(count),embryo);
        cousin_metrics = [cousin_metrics;out1,out2,out3,out4,out5,out6,out7];
    end
    normed_cousin_metrics = zscore(cousin_metrics);
    
    %end result:zscore of the duration of cell i's cell cycle compared to
    %ancestors'/cousins'
    similarity_to_ancestors = normed_ancestor_metrics(end,3); %hard-coded to grab duration of cell cycle column
    similarity_to_cousins = normed_cousin_metrics(end,3);
    

    
    
    
    
    
    
    
    
    
    
    %%
    %Cell density
    
    %hyperparamters
    %radius of sphere--relates to axis scale; implement conversion later
    r = 20; %microns
    r = r/.16;%convert using anisotropy
    sigma = r/3;
    
    
    cell_position = embryo(t).finalpoints(i,:);
    in_sphere_inds = [];
    for j = 1:n
        if norm(embryo(t).finalpoints(j,:)-cell_position) <= r
            in_sphere_inds = [in_sphere_inds;j];
        end
    end
    
    %make cell i origin for points in the sphere
    in_sphere_points = embryo(t).finalpoints - cell_position;
    
    first_deriv_x = 0;
    first_deriv_y = 0;
    first_deriv_z = 0;
    laplacian = 0;
    
    
    for j = 1:size(in_sphere_points,1)
        xx = in_sphere_points(j,1);
        yy = in_sphere_points(j,2);
        zz = in_sphere_points(j,3);
        
        %CHECK SIGN CORRECTNESS/CONSISTENCY W DESIRED CONVENTION IN
        %GAUSSIANS FUNCTION
        first_deriv_x = first_deriv_x + gaussians(xx,yy,zz,sigma,"first_x");
        first_deriv_y = first_deriv_y + gaussians(xx,yy,zz,sigma,"first_y");
        first_deriv_z = first_deriv_z + gaussians(xx,yy,zz,sigma,"first_z");
        laplacian = laplacian + gaussians(xx,yy,zz,sigma,"laplacian");
    end
    %now have the weighted anisotropies of cell density in x,y,z directions
    %and laplacian
    
    %abstract some of the above into separate functions
    
    
    
    
    
    
    
    
    
    
    %%
   %Local normal
    
    %hyperparameters are knorm and k
    k_residual = 5;
    %knorm = 5; (implemented above outside of loop)
    
    %local normal vector
    local_normal_vector = norms(i,:);
    
    %spatial robustness of normal
    % distance to k nns
    [~,indices_of_residual_nns] = pdist2(points,points(i,:),'euclidean','Smallest',k_residual+1);
    indices_of_residual_nns = indices_of_residual_nns(2:end);
    
    vecs_to_nns = embryo(t).finalpoints(indices_of_residual_nns);
    vecs_to_nns = vecs_to_nns - embryo(t).finalpoints(i,:);
    
    residual = 0;
    for j = 1:size(vecs_to_nns,1)
        residual = residual + abs(dot(vecs_to_nns(j,:),local_normal_vector));
    end
    
    
    
    
    
    
    
    
    
    
    
    
    %%
    %Inter-cell distance
    
    %hyperparameters
    k = 10; %number of neighbors %%%CURRENTLY HAS TO BE SAME AS k_residual; CHANGE!!!!!
    
    %distance to k nns
    [distances_to_k_nns,indices_of_k_nns] = pdist2(points,points(i,:),'euclidean','Smallest',k+1);
    distances_to_k_nns = distances_to_k_nns(2:end,:);
    %make horizontal (for concatenation below)
    distances_to_k_nns = distances_to_k_nns';
    
    indices_of_k_nns = indices_of_k_nns(2:end);
    
    
    %signed distances
    signed_distances = [];
    for idx = 1:size(vecs_to_nns,1)
        signed_distances = [signed_distances,vecs_to_nns(idx,:)];
    end

    
    
    
    
    
    
    
    
    
    
   %% 
   %Migration paths
    
    %hyperparameter mins (each feature could have its own time window)
    mins = 10;%minutes
    %convert
    mins = ceil(mins*305/470); %305 timesteps before twitching/470 mins (on average)
    
    
    %vec connecting cell to where it was x mins ago
    [past_time,past_ind,~,~] = backTraverse(t,i,'timesteps',mins,embryo);
    vec1 = embryo(past_time).finalpoints(past_ind,:) - embryo(t).finalpoints(i,:);
    
    %vec connecting cell to where it will be x mins in future (or end of dataset)
    [future_pts] = forwardTraverse(t,[i],'timesteps',t+mins,embryo);
    
    positions = [];
    
    for j = 1:length(future_pts)
        positions = [positions;embryo(min(t+mins,length(embryo)-1)).finalpoints(future_pts(j),:)];
    end
    %in case cell/lineage dies
    if isempty(positions)
        vec2 = [0,0,0];
    else
        positions = mean(positions,1);
        vec2 = positions - embryo(t).finalpoints(i,:);
    end
        

    
    %vec connecting cell to site of last division
    [past_division_time,past_division_ind,~,~] = backTraverse(t,i,'generations',1,embryo);
    vec3 = embryo(past_division_time).finalpoints(past_division_ind,:) - embryo(t).finalpoints(i,:);
    
    %vec connecting cell to site of next division (0 if cell dies)
    timept = t;
    indpt = i;
    vec4 = [];
    if embryo(timept).suc(indpt,1) == -1
        vec4 = [0,0,0];
    else
        while embryo(timept).suc(indpt,2) == -1
            indpt = embryo(timept).suc(indpt,1);
            timept = timept + 1;
            %if at end of embryo abort
            if timept >= length(embryo)-1
                vec4 = [0,0,0];
                break
            end
            %abort if cell dies (at current or next time)
            if embryo(timept).suc(indpt,1) == -1
                vec4 = [0,0,0];
                break
            end
        end
    end
    
    if isempty(vec4)
        vec4 = embryo(timept).finalpoints(indpt,:) - embryo(t).finalpoints(i,:);
    end
    
    
    %vec connecting site of last division to site of next division
    %uses previously defined quantities
    
    %handle having reached the end of embryo
    if timept > length(embryo)-1
        vec5 = [0,0,0];
    else
        vec5 = embryo(timept).finalpoints(indpt,:) - embryo(past_division_time).finalpoints(past_division_ind,:);
    end
    
    
    
    
    
    
    
    
    
    
    
    
    %%
    %congruency with neighbors' motion
    
    %hyperparameters
    %minutes_back
    %k_congruency
    %k_congruency_present_time
    
    %go back in time
    minutes_back = 10;%mins
    %convert to timesteps
    minutes_back = ceil(minutes_back*305/470);
    
    [previous_time,previous_ind,~,~] = backTraverse(t,i,'timesteps',minutes_back,embryo);
    
    %compute nns
    k_congruency = 40;
    [~,indices_of_congruency_nns] = pdist2(embryo(previous_time).finalpoints,embryo(previous_time).finalpoints(previous_ind,:),'euclidean','Smallest',k_congruency+1);
    indices_of_congruency_nns = indices_of_congruency_nns(2:end);
    
    %advance cell i's ancestor and its nns to present time
    descendants_of_cells_ancestor = forwardTraverse(previous_time,previous_ind,'timesteps',t,embryo);
    descendants_of_past_neighbors = forwardTraverse(previous_time,indices_of_congruency_nns,'timesteps',t,embryo);
    
    
    if isempty(descendants_of_cells_ancestor) || isempty(descendants_of_past_neighbors)
        present_nns_that_are_descendants_of_past_nns = 0;
    else
        %calculate nns of current descendants of ancestor cell
        k_congruency_present_time = 40;

        present_nns = [];
        for j = 1:length(descendants_of_cells_ancestor)
            [~,present_inds] = pdist2(embryo(t).finalpoints,embryo(t).finalpoints(descendants_of_cells_ancestor(j),:),'euclidean','Smallest',k_congruency_present_time+1);
            present_inds = present_inds(2:end);

            present_nns = union(present_nns,present_inds);
        end

        %calculate proportion of present nns that are descendants of ancestor
        proportion_of_present_nns_that_are_descendants_of_past_nns = sum(ismember(present_nns,descendants_of_past_neighbors))/length(present_nns);
    end
    
    
    
    
    
    
    
    
    
    
    
    %%
    %{
    %this is now obsolete
    %main path
    %strategy 1--fuzzy graph
    
    %hyper parameter- amount of time back- using mins from above, but could
    %be different
    
    %fuzzy radius 
    fuzzy_radius = 10; %this will need to be converted to coordinate space of emb dataset
    
    
    tpt = t;
    it = i;

    counter = mins;
    
    %generate list of positions
    path_positions = [embryo(tpt).finalpoints(ipt,:)];
    while counter > 0
        %walk back
        ind_pt = stepback(time_pt,ind_pt,embryo);
        time_pt = time_pt-1;
        %log position of cell at earlier time
        path_positions = [path_positions;embryo(tpt).finalpoints(ipt,:)];
        counter = counter-1;
    end
    
    %make adjacency matrix
    adj = zeros(size(path_positions,1));
    for j = 1:size(path_positions,1)-1
        adj(j,j+1) = norm(path_positions(j,:)-path_positions(j+1,:));
        adj(j+1,j) = norm(path_positions(j,:)-path_positions(j+1,:));
    end
    
    %add connections by spatial proximity
    for j = 1:length(path_positions)
        for h = 1:length(path_positions)
            d =  norm(path_positions(j,:) - path_positions(h,:));
            if d <= fuzzy_radius
                adj(j,h) = d;
                adj(h,j) = d;
            end
        end
    end
    
    fuzzy_graph = graph(adj);
    shortest_path = shortestpath(fuzzy_graph,1,size(path_positions,1));
    
    %ratio of length of shortestpath,displacement
    
    shortest_path_length = 0;
    for j = 1:length(shortest_path-1)
        shortest_path_length = shortest_path_length + norm(path_positions(j,:),path_positions(j+1));
    end
    
    directness_ratio = norm(path_positions(1,:)-path_positions(size(path_positions,1),:))/shortest_path_length;
    %}
    
    
    
    
    
    
    
    
    
    
    
    %%
    %nuclear size
    raw_nuc_size = embryo(t).celldata(i,7); %just happens to be how it's loaded from the disk
    
    r_nuc = 15; %microns
    %convert to unites using anisotropy
    r_nuc = r_nuc/.16;
    
    nuc_neighbors = [];
    for j = 1:length(embryo(t).finalpoints)
        if norm(embryo(t).finalpoints(i,:) - embryo(t).finalpoints(j,:)) < r_nuc
            nuc_neighbors = [nuc_neighbors;j];
        end
    end
    
    filtered_nuc_size = mean(embryo(t).celldata(nuc_neighbors,7));%includes cell i
    
    %{
    k_nuc = 10;
    [~,indices_of_nuc_nns] = pdist2(embryo(t).finalpoints,embryo(t).finalpoints(i,:),'euclidean','Smallest',k_nuc+1);
    
    filtered_nuc_size = mean(embryo(t).celldata(indices_of_nuc_nns,7));
    %}
    
    

    
    
    
    
    %%
    %developmental stage
    %number of cells in embryo
    number_of_cell_in_embryo = n;
    
    %developmental time (relative to 4 cell stage)***adjust to account for
    %temp/UV light--involves framerate (part of dataset)
    developmental_time = t;
    %accounting for 4cell stage time is done in normalization (normalize3)
    
    
    
    
    
    
 
    %%

    %%% Feature Vector %%%
    
    
    featurevector = [featurevector;beginning_of_current_cell_cycle,end_of_cell_current_cycle,end_of_cell_current_cycle,time_until_next_division,time_since_last_division,...
                                    fraction_of_current_cycle_elapsed,similarity_to_ancestors,similarity_to_cousins,first_deriv_x,first_deriv_y,first_deriv_z,laplacian,...
                                    local_normal_vector,residual,distances_to_k_nns,signed_distances,vec1,vec2,vec3,vec4,vec5,proportion_of_present_nns_that_are_descendants_of_past_nns,...
                                    filtered_nuc_size,number_of_cell_in_embryo,developmental_time];
   
    
    
    
    
    
    %{ 
    %vertical formatting doesn't work for some reason, so it's here for convenience                                
    featurevector = [featurevector;beginning_of_current_cell_cycle,
                                    end_of_cell_current_cycle,
                                    duration_of_current_cell_cycle,
                                    time_until_next_division,
                                    time_since_last_division,
                                    fraction_of_current_cycle_elapsed,
                                    similarity_to_ancestors,
                                    similarity_to_cousins,
                                    first_deriv_x,
                                    first_deriv_y,
                                    first_deriv_z,
                                    laplacian,
                                    local_normal_vector,
                                    residual,
                                    distances_to_k_nns,
                                    signed_distances,
                                    vec1,
                                    vec2,
                                    vec3,
                                    vec4,
                                    vec5,
                                    proportion_of_present_nns_that_are_descendants_of_past_nns,
                                    filtered_nuc_size,
                                    number_of_cell_in_embryo,
                                    developmental_time];
    %}
    
    
    
    
    
    
    
    
    
    
    %all 
    %featurevector = [featurevector;distfeat',normprojfeat',normscalproj',perpprojfeat',perpscalproj',rel_position,num_cells,distcos'];
    %lacking rel_position
    %featurevector = [featurevector;distfeat',normprojfeat',normscalproj',perpprojfeat',perpscalproj',num_cells,distcos'];
    
    %only local
    %featurevector = [featurevector;distfeat',normprojfeat',normscalproj',perpprojfeat',perpscalproj',distcos'];
    
    %for fish
    %featurevector = [featurevector;distfeat',normprojfeat',normscalproj',perpprojfeat',perpscalproj'];
    
    
    %looking for branching
    %featurevector = [featurevector;distfeat',rel_position,num_cells];
    
    
    
    
    
    %used to evaluate effect of alignment--has no global features
    %featurevector = [featurevector;distfeat',normprojfeat',normscalproj',perpprojfeat',perpscalproj',num_cells];
    
    
    %featurevector = [featurevector;motion];
    
    
    
    %featurevector = [featurevector;distfeat',normprojfeat',time];
    %featurevector = [featurevector;normprojfeat'];
    
    %5 best local features
    %featurevector = [featurevector;distfeat',normprojfeat',normprojfeat',perpprojfeat',perpscalproj'];
    
    
    %{
    %%% Number of elements of each feature %%% *** automate this
    numElms = [k,3*k,k,3*k,k,3,1,3];
    featTypes = ['local','local3d','local','local3d','local','3d','singlevalue'];
    %TODO--make it so the normalization goes according to this, automate
    %these additional outputs
    %}
    
    %note: all vectors passed into this vector must be horizontal
    %note: each row in featurevector is one particular cell at one particular time
    
    
    %%% Helping functions %%%
    
    
    
    
    
    

end

end
