%NOW OUTDATED!!!!!!!!!

function [featurevector] = computeFeaturesTimepoint3(embryo,t,feats_and_norms,param_values)
%compute feature vector on each point in individual point cloud;


%defining parameters
%{
radius_bounding_sphere (density gradient and nuclear size)
k for local normal calc and residual
k for dists/vecs to nns
degree of cousins
amount of time for vecs to positions in past/future
amount of time for frac of nns that remain nns
number of nns that remain nns
%}

%would be better if each cell entry were an array of the hyperparams and
%then the next entry were an array of values (how to work in multiple values?)
%make finding whether points are within bounding sphere more efficient



%make these actually appropriate
default_param_values = {"radius_boundingsphere",            10;
                        "k_intercelldistance",              25;
                        "degree_cousins",                   1;
                        "timewindow_migrationpaths",        10;
                        "timewindow_migrationcongruency",   30;
                        "k_migrationcongruency",            10};
 
 
%set hyperparameters
if isempty(find([param_values{:,1}]=="radius_boundingsphere"))
    radius_boundingsphere = default_param_values{([default_param_values{:,1}]=="radius_boundingsphere"),2};
else
    radius_boundingsphere = param_values{([param_values{:,1}]=="radius_boundingsphere"),2};
end

if isempty(find([param_values{:,1}]=="k_intercelldistance"))
    k_intercelldistance = default_param_values{([default_param_values{:,1}]=="k_intercelldistance"),2};
else
    k_intercelldistance = param_values{([param_values{:,1}]=="k_intercelldistance"),2};
end

if isempty(find([param_values{:,1}]=="degree_cousins"))
    degree_cousins = default_param_values{([default_param_values{:,1}]=="degree_cousins"),2};
else
    degree_cousins = param_values{([param_values{:,1}]=="degree_cousins"),2};
end

if isempty(find([param_values{:,1}]=="timewindow_migrationpaths"))
    timewindow_migrationpaths = default_param_values{([default_param_values{:,1}]=="timewindow_migrationpaths"),2};
else
    timewindow_migrationpaths = param_values{([param_values{:,1}]=="timewindow_migrationpaths"),2};
end

if isempty(find([param_values{:,1}]=="timewindow_migrationcongruency"))
    timewindow_migrationcoherency = default_param_values{([default_param_values{:,1}]=="timewindow_migrationcongruency"),2};
else
    timewindow_migrationcoherency = param_values{([param_values{:,1}]=="timewindow_migrationcongruency"),2};
end

if isempty(find([param_values{:,1}]=="k_migrationcongruency"))
    k_migrationcongruency = default_param_values{([default_param_values{:,1}]=="k_migrationcongruency"),2};
else
    k_migrationcongruency = param_values{([param_values{:,1}]=="k_migrationcongruency"),2};
end

%{
radius_boundingsphere = default_param_values{([default_param_values{:,1}=="radius_boundingsphere"]),2};
k_normal = default_param_values{([default_param_values{:,1}=="k_normal"]),2};
k_intercelldistance = default_param_values{([default_param_values{:,1}=="k_intercelldistance"]),2};
degree_cousins = default_param_values{([default_param_values{:,1}=="degree_cousins"]),2};
timewindow_migrationpaths = default_param_values{([default_param_values{:,1}=="timewindow_migrationpaths"]),2};
timewindow_migrationcoherency = default_param_values{([default_param_values{:,1}=="timewindow_migrationcoherency"]),2};
k_migrationcoherency = default_param_values{([default_param_values{:,1}=="k_migrationcoherency"]),2};
%}


%features

%this list contains every feature that can possibly be added
ordered_feats = ["beginning_of_current_cell_cycle"; 
              "end_of_cell_current_cycle";       
              "duration_of_current_cell_cycle";  
              "time_until_next_division";        
              "time_since_last_division";        
              "fraction_of_current_cycle_elapsed";
              "similarity_to_ancestors";          
              "similarity_to_cousins";            
              "first_deriv_x";                    
              "first_deriv_y";                    
              "first_deriv_z";                    
              "laplacian";                        
              "local_normal_vector";              
              "residual";                         
              "distances_to_k_nns"                  
              "signed_distances"                  
              "vec1";                             
              "vec2";                             
              "vec3";                            
              "vec4";                             
              "vec5";                             
              "proportion_of_present_nns_that_are_descendants_of_past_nns"; 
              "filtered_nuc_size";                
              "number_of_cell_in_embryo";         
              "developmental_time"];    
%populate indicators based on user input
indicators = zeros(size(ordered_feats));
for i = 1:length(indicators)
    indicators(i) = double(~isempty(find(feats_and_norms(:,1)==ordered_feats(i))));
end


%ordered_featuresizes







%aux for the computation

%points from embryo--non-temporal features work within one time step
points = embryo(t).finalpoints;

%num points
n = size(points,1);

%normals, calculated with k_norm nns
[~,normals,residuals] = cohereNormals(points,radius_boundingsphere/.16,3); % 
%NOTE last argument of cohere normals is c = interconnectivity of graph


featurevector = [];
featuresizes = [];




%features

for i=1:n
    
    cell_featurevector = [];
    
    %%
    %Cell Cycle
    
    %if any of the cell cycle metrics is requested, compute
    if sum(indicators(1:6)>0)

        %cell cycle metrics of current cell
        [beginning_of_current_cell_cycle,end_of_cell_current_cycle,duration_of_current_cell_cycle,time_until_next_division,time_since_last_division,...
            fraction_of_current_cycle_elapsed,time_of_next_division] = cellCycleMetrics(t,i,embryo);
        %leave out time of next division; redundant
        
        %add desired quantities to cell_featurevector
        if indicators(ordered_feats=="beginning_of_current_cell_cycle")
            cell_featurevector = [cell_featurevector,beginning_of_current_cell_cycle];
        end
        if indicators(ordered_feats=="end_of_cell_current_cycle")
            cell_featurevector = [cell_featurevector,end_of_cell_current_cycle];
        end
        if indicators(ordered_feats=="duration_of_current_cell_cycle")
            cell_featurevector = [cell_featurevector,duration_of_current_cell_cycle];
        end
        if indicators(ordered_feats=="time_until_next_division")
            cell_featurevector = [cell_featurevector,time_until_next_division];
        end
        if indicators(ordered_feats=="time_since_last_division")
            cell_featurevector = [cell_featurevector,time_since_last_division];
        end
        if indicators(ordered_feats=="fraction_of_current_cycle_elapsed")
            cell_featurevector = [cell_featurevector,fraction_of_current_cycle_elapsed];
        end

    end

    
    %%
    if sum(indicators(7:8)>0)
        %Comparison with Cousins

        %hyperparameter
        %number of generations back to go (in this case both to find ancestors and cousins)
        g = degree_cousins;

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
        %manually extract only most recent ancestor
        ancestor_metrics = [ancestor_metrics(1,:);ancestor_metrics(end,:)];

       
        %cousins
        cousin_indices = [cousin_indices,i]; %i in last row
        %remove any dead cells that may have entered the cousins
        cousin_indices(find(cousin_indices==-1)) = [];

        cousin_metrics = [];
        for count = 1:length(cousin_indices)
            [out1,out2,out3,out4,out5,out6,out7] = cellCycleMetrics(t,cousin_indices(count),embryo);
            cousin_metrics = [cousin_metrics;out1,out2,out3,out4,out5,out6,out7];
        end
        
        %{
        %debugging
        i
        'cell i'
        embryo(t).cellnames(i)
        'cousins'
        embryo(t).cellnames(cousin_indices)
        cousin_metrics(:,3)
        %}
        
        if size(cousin_metrics,1) == 1
            similarity_to_cousins = 1;
        else %assuming it's never zero
            similarity_to_cousins = cousin_metrics(end,3)/(mean(cousin_metrics(1:end-1,3))); %3 is location of cell cycle
            %ratio between cell i's cell cycle duration and cousins' cell
            %cycle durations
        end
        
        if size(ancestor_metrics,1) == 1
            similarity_to_ancestors = 1;
        else %assuming it's never zero
            similarity_to_ancestors = ancestor_metrics(end,3)/(mean(ancestor_metrics(1:end-1,3))); %3 is location of cell cycle
            %ratio between cell i's cell cycle duration and ancestors' cell
            %cycle durations
        end

        
        %add desired quantities to cell_featurevector
        if indicators(ordered_feats=="similarity_to_ancestors")
            cell_featurevector = [cell_featurevector,similarity_to_ancestors];
        end
        if indicators(ordered_feats=="similarity_to_cousins")
            cell_featurevector = [cell_featurevector,similarity_to_cousins];
        end
        
    end
    
    %%
    %Cell density
    
    if sum(indicators(9:12)>0)
    
        %hyperparamters
        %radius of sphere--relates to axis scale; implement conversion later
        r = radius_boundingsphere; %microns
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

            first_deriv_x = first_deriv_x + gaussians(xx,yy,zz,sigma,"first_x");
            first_deriv_y = first_deriv_y + gaussians(xx,yy,zz,sigma,"first_y");
            first_deriv_z = first_deriv_z + gaussians(xx,yy,zz,sigma,"first_z");
            laplacian = laplacian + gaussians(xx,yy,zz,sigma,"laplacian");
        end
        %now have the weighted anisotropies of cell density in x,y,z directions
        %and laplacian
        
        %add desired quantities to cell_featurevector
        if indicators(ordered_feats=="first_deriv_x")
            cell_featurevector = [cell_featurevector,first_deriv_x];
        end
        if indicators(ordered_feats=="first_deriv_y")
            cell_featurevector = [cell_featurevector,first_deriv_y];
        end
        if indicators(ordered_feats=="first_deriv_z")
            cell_featurevector = [cell_featurevector,first_deriv_z];
        end
        if indicators(ordered_feats=="laplacian")
            cell_featurevector = [cell_featurevector,laplacian];
        end
    end
    
    
    
    %%
    %Local normal
    
   
    if sum(indicators(13:14)>0)

        %local normal vector
        local_normal_vector = normals(i,:);
        
        %residual
        residual = residuals(i);
        
        %add desired quantities to cell_featurevector
        if indicators(ordered_feats=="local_normal_vector")
            cell_featurevector = [cell_featurevector,local_normal_vector];
        end
        if indicators(ordered_feats=="residual")
            cell_featurevector = [cell_featurevector,residual];
        end

    end
    
    

    
    
    %%
    %Inter-cell distance
    
    if sum(indicators(15:16)>0)
    
        %hyperparameters
        k = k_intercelldistance;

        %distance to k nns
        [distances_to_k_nns,indices_of_k_nns] = pdist2(embryo(t).finalpoints,embryo(t).finalpoints(i,:),'euclidean','Smallest',k+1);
        distances_to_k_nns = distances_to_k_nns(2:end,:);
        %make horizontal (for concatenation below)
        distances_to_k_nns = distances_to_k_nns';

        indices_of_k_nns = indices_of_k_nns(2:end);
        vecs_to_nns = embryo(t).finalpoints(indices_of_k_nns,:) - embryo(t).finalpoints(i,:);

        %signed distances
        signed_distances = [];
        for idx = 1:size(vecs_to_nns,1)
            signed_distances = [signed_distances,vecs_to_nns(idx,:)];
        end
        
        %add desired quantities to cell_featurevector
        if indicators(ordered_feats=="distances_to_k_nns")
            cell_featurevector = [cell_featurevector,distances_to_k_nns];
        end
        if indicators(ordered_feats=="signed_distances")
            cell_featurevector = [cell_featurevector,signed_distances];
        end

    end
    
    

    
    
    
    %% 
    %Migration paths **RENAME vec1...evc5 including the strings
    
    if sum(indicators(17:21)>0)
   
        %hyperparameter mins (each feature could have its own time window)
        mins = timewindow_migrationpaths;%minutes
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
        
        %add desired quantities to cell_featurevector
        if indicators(ordered_feats=="vec1")
            cell_featurevector = [cell_featurevector,vec1];
        end
        if indicators(ordered_feats=="vec2")
            cell_featurevector = [cell_featurevector,vec2];
        end
        if indicators(ordered_feats=="vec3")
            cell_featurevector = [cell_featurevector,vec3];
        end
        if indicators(ordered_feats=="vec4")
            cell_featurevector = [cell_featurevector,vec4];
        end
        if indicators(ordered_feats=="vec5")
            cell_featurevector = [cell_featurevector,vec5];
        end
        
    end
    
    
    
    
    %%
    %congruency with neighbors' motion
    
    if indicators(22)

        %hyperparameters
        %minutes_back
        %k_congruency
        %k_congruency_present_time

        %go back in time
        minutes_back = timewindow_migrationcoherency;%mins
        %convert to timesteps
        minutes_back = ceil(minutes_back*305/470);

        [previous_time,previous_ind,~,~] = backTraverse(t,i,'timesteps',minutes_back,embryo);

        %compute nns
        k_congruency = k_migrationcongruency;
        [~,indices_of_congruency_nns] = pdist2(embryo(previous_time).finalpoints,embryo(previous_time).finalpoints(previous_ind,:),'euclidean','Smallest',k_congruency+1);
        indices_of_congruency_nns = indices_of_congruency_nns(2:end);

        %advance cell i's ancestor and its nns to present time
        descendants_of_cells_ancestor = forwardTraverse(previous_time,previous_ind,'timesteps',t,embryo);
        descendants_of_past_neighbors = forwardTraverse(previous_time,indices_of_congruency_nns,'timesteps',t,embryo);


        if isempty(descendants_of_cells_ancestor) || isempty(descendants_of_past_neighbors)
            present_nns_that_are_descendants_of_past_nns = 0;
        else
            %calculate nns of current descendants of ancestor cell
            k_congruency_present_time = k_migrationcongruency;

            present_nns = [];
            for j = 1:length(descendants_of_cells_ancestor)
                [~,present_inds] = pdist2(embryo(t).finalpoints,embryo(t).finalpoints(descendants_of_cells_ancestor(j),:),'euclidean','Smallest',k_congruency_present_time+1);
                present_inds = present_inds(2:end);

                present_nns = union(present_nns,present_inds);
            end

            %calculate proportion of present nns that are descendants of ancestor
            proportion_of_present_nns_that_are_descendants_of_past_nns = sum(ismember(present_nns,descendants_of_past_neighbors))/length(present_nns);
        end
        
        %add desired quantity to cell_featurevector (only feature for this category, so no need to check)
        cell_featurevector = [cell_featurevector,proportion_of_present_nns_that_are_descendants_of_past_nns];
        
    end
    
    
    
    
    
    
    %%
    %nuclear size
    
    if indicators(23)
    
        raw_nuc_size = embryo(t).celldata(i,7); %just happens to be how it's loaded from the disk

        r_nuc = radius_boundingsphere; %microns
        %convert to unites using anisotropy
        r_nuc = r_nuc/.16;

        nuc_neighbors = [];
        for j = 1:length(embryo(t).finalpoints)
            if sum((embryo(t).finalpoints(i,:) - embryo(t).finalpoints(j,:)).^2) < r_nuc^2
                nuc_neighbors = [nuc_neighbors;j];
            end
        end

        filtered_nuc_size = mean(embryo(t).celldata(nuc_neighbors,7));%includes cell i
        
        cell_featurevector = [cell_featurevector,filtered_nuc_size];
    
    end
    

    
    
    
    
    %%
    %developmental stage
    
    if indicators(24)
    
        %number of cells in embryo
        number_of_cell_in_embryo = n;
    
        %add desired quantity to cell_featurevector (only feature for this category, so no need to check)
        cell_featurevector = [cell_featurevector,number_of_cell_in_embryo];
        
    end
    
    if indicators(24)
        
        %developmental time (relative to 4 cell stage)***adjust to account for
        %temp/UV light--involves framerate (part of dataset)
        developmental_time = t;
        %accounting for 4cell stage time is done in normalization (normalize3)
    
        %add desired quantity to cell_featurevector (only feature for this category, so no need to check)
        cell_featurevector = [cell_featurevector,developmental_time];
        
    end
    
    
    
 
    %%

    %%% Feature Vector %%%
    
    
    featurevector = [featurevector;cell_featurevector];
   
    
    
    
end

end
