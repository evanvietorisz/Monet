function [featurevector,normed_featurevector,indices_list] = computeFeatures4(embryo,times,feats,param_values,normalize_bool,input_xy_res)
%{
%COMPUTEFEATURES4
    Compute and normalize spatiotemporal feature vector for each cell at each time in an embyro struct


ARGUMENTS
_________
embryo (struct array)
    the embryo struct from which to compute features

times (array)
    the times (frames) from embryo struct at which cells will have their
    feature representations computed

feats (string array)
    contains the names of the features you wish to compute; can be in any
    order

param_values (cell array)
    contains the names (column 1) and values (column 2) of the
    hyperparameters required to compute features; can be in any order

normalize_bool (bool)
    whether or not to normalize

input_xy_res (scalar)
    xy resolution of the embryo


OUTPUTS
_______
featurevector (array)
    matrix in which rows are feature represnetations of each indiudal cell,
    unnormalized

normed_featurevector (array) 
    same as featurevector but normalized according to the scheme below.
    features as listed in ordered_feats correspond to the normalizations in
    ordered_norms

indices_list (array of shape # features x 2)
    array whose columns contain the starting and ending indices of each
    feature, respectively. ex. if feature 3 goes from columns 6 to 8 in the
    featurevector, then indices_list(3,:) = [6,8]
%}



default_param_values = {"radius_boundingsphere_celldensity",20;
                        "radius_boundingsphere_residual",   20;
                        "radius_boundingsphere_nucsize",    20;
                        "k_intercelldistance",              50;
                        "degree_cousins",                   1;
                        "timewindow_migrationpaths",        .125;
                        "timewindow_migrationcongruency",   .125;
                        "k_migrationcongruency",            50};
                    
                    
params_and_corresponding_feats = {"radius_boundingsphere_celldensity",["first_deriv_x";"first_deriv_y";"first_deriv_z";"laplacian"];
                                  "radius_boundingsphere_residual",   ["local_normal_vector";"residual"];
                                  "radius_boundingsphere_nucsize",    ["filtered_nuc_size"];
                                  "k_intercelldistance",              ["distances_to_k_nns";"signed_distances"];
                                  "degree_cousins",                   ["similarity_to_ancestors_cell_cycle";"similarity_to_cousins_cell_cycle"]
                                  "timewindow_migrationpaths",        ["vec_to_position_x_mins_ago","vec_to_position_x_mins_in_future"];
                                  "timewindow_migrationcongruency",   ["proportion_of_present_nns_that_are_descendants_of_past_nns"];
                                  "k_migrationcongruency",            ["proportion_of_present_nns_that_are_descendants_of_past_nns"]};
                    
         
                    
 
%set hyperparameters

if isempty(find([param_values{:,1}]=="radius_boundingsphere_celldensity"))
    radius_boundingsphere_celldensity = default_param_values{([default_param_values{:,1}]=="radius_boundingsphere_celldensity"),2};
else
    radius_boundingsphere_celldensity = param_values{([param_values{:,1}]=="radius_boundingsphere_celldensity"),2};
end

if isempty(find([param_values{:,1}]=="radius_boundingsphere_residual"))
    radius_boundingsphere_residual = default_param_values{([default_param_values{:,1}]=="radius_boundingsphere_residual"),2};
else
    radius_boundingsphere_residual = param_values{([param_values{:,1}]=="radius_boundingsphere_residual"),2};
end

if isempty(find([param_values{:,1}]=="radius_boundingsphere_nucsize"))
    radius_boundingsphere_nucsize = default_param_values{([default_param_values{:,1}]=="radius_boundingsphere_nucsize"),2};
else
    radius_boundingsphere_nucsize = param_values{([param_values{:,1}]=="radius_boundingsphere_nucsize"),2};
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
    timewindow_migrationcongruency = default_param_values{([default_param_values{:,1}]=="timewindow_migrationcongruency"),2};
else
    timewindow_migrationcongruency = param_values{([param_values{:,1}]=="timewindow_migrationcongruency"),2};
end

if isempty(find([param_values{:,1}]=="k_migrationcongruency"))
    k_migrationcongruency = default_param_values{([default_param_values{:,1}]=="k_migrationcongruency"),2};
else
    k_migrationcongruency = param_values{([param_values{:,1}]=="k_migrationcongruency"),2};
end



%features

%this list contains every feature that can possibly be added
ordered_feats = ["beginning_of_current_cell_cycle"; 
              "end_of_current_cell_cycle";       
              "duration_of_current_cell_cycle";  
              "time_until_next_division";        
              "time_since_last_division";        
              "fraction_of_current_cycle_elapsed";
              "similarity_to_ancestors_cell_cycle";          
              "similarity_to_cousins_cell_cycle";            
              "first_deriv_x";                    
              "first_deriv_y";                    
              "first_deriv_z";                    
              "laplacian";                        
              "local_normal_vector";              
              "residual";                         
              "distances_to_k_nns"                  
              "signed_distances"                  
              "vec_to_position_x_mins_ago";                             
              "vec_to_position_x_mins_in_future";                             
              "vec_to_site_of_last_division";                            
              "vec_to_site_of_next_division";                             
              "displacement_vec_over_cell_cycle";                             
              "proportion_of_present_nns_that_are_descendants_of_past_nns"; 
              "filtered_nuc_size";                
              "number_of_cell_in_embryo";         
              "developmental_time"];    

%normalizations corresponding to features above
ordered_norms = ["scale_4cell_to_end"; 
              "scale_4cell_to_end";       
              "scale_4cell_to_end";  
              "scale_4cell_to_end";        
              "scale_4cell_to_end";        
              "nothing";
              "nothing";          
              "nothing";            
              "column_zscore";                    
              "column_zscore";                    
              "column_zscore";                    
              "column_zscore";                        
              "divide_avg_mag";              
              "column_zscore";                         
              "all_zscore"                  
              "divide_avg_mag"                  
              "divide_avg_mag";                             
              "divide_avg_mag";                             
              "divide_avg_mag";                            
              "divide_avg_mag";                             
              "divide_avg_mag";                             
              "nothing"; 
              "column_zscore";                
              "scale_0_to_max_num_cells";         
              "scale_4cell_to_end"];      
      
          
%for computing what indices a feature corresponds to
ordered_featuresizes = [1;
                        1;
                        1;
                        1;
                        1;
                        1;
                        1;
                        1;
                        1;
                        1;
                        1;
                        1;
                        3;
                        1;
                        k_intercelldistance;
                        3*k_intercelldistance;
                        3;
                        3;
                        3;
                        3;
                        3;
                        1;
                        1;
                        1;
                        1];
            
                    
%populate indicators based on user input
indicators = zeros(size(ordered_feats));
for i = 1:length(indicators)
    indicators(i) = double(~isempty(find(feats(:,1)==ordered_feats(i))));
end


%make indices_list
kept_sizes = ordered_featuresizes(indicators==1);%actually requested by user
kept_features = ordered_feats(indicators==1);

kept_indices_list = zeros(size(feats,1),2);
c = 1;
for i = 1:size(feats,1)
    kept_indices_list(i,1) = c;
    kept_indices_list(i,2) = c + kept_sizes(i)-1;
    c = c + kept_sizes(i);
end

indices_list = zeros(size(kept_indices_list));
for i = 1:size(kept_indices_list,1)
    indices_list(i,:) = kept_indices_list(kept_features==feats(i,1),:);
end

%convert hyperparameters into native units

%time windows
%   parameter value t is fraction of time from 4 cell to twitch ~~ 400 mins
%   eg. 50 mins = .125*400, so 50 mins should be input as .125
%   #tpts = ceil(t*#tpts_4cell_to_twitch)
%   twitch may not be measurable in dataset, so approx using 350cell tpt
%   #ftpts_4cell_to_350cell = p*(#tpts_4cell_to_350cell)
%   p = approximate ratio of time from 4cell to twitch over 4cell to 350
%   cell = (475-75)/(345-75)  (from https://www.wormatlas.org/hermaphrodite/introduction/mainframe.html)

p = (475-75)/(345-75);

%obtain first time at which emb is 4 cells
fourcelltime = 1;
while size(embryo(fourcelltime).finalpoints,1) < 4
    fourcelltime = fourcelltime + 1;
end

%obtain first time at which emb is 350 cells
threefiftycelltime = 1;
while size(embryo(threefiftycelltime).finalpoints,1) < 350
    threefiftycelltime = threefiftycelltime + 1;
end

%approximate twitch tpt
%   approx_twitch_time - fourcelltime = p*(threefiftycelltime-fourcelltime)
%   p = approximate ratio of time from 4cell to twitch over 4cell to 350
%   cell = (475-75)/(345-75)  (from https://www.wormatlas.org/hermaphrodite/introduction/mainframe.html)

approx_twitch_time = p*(threefiftycelltime - fourcelltime) + fourcelltime;

%convert raw parameter values to corresponding numbers of tpts

timewindow_migrationpaths = ceil(timewindow_migrationpaths * (approx_twitch_time-fourcelltime));
timewindow_migrationcongruency = ceil(timewindow_migrationcongruency * (approx_twitch_time-fourcelltime));




%do conversion for input_xy_res, take out of below
radius_boundingsphere_celldensity = radius_boundingsphere_celldensity/input_xy_res;
radius_boundingsphere_residual = radius_boundingsphere_residual/input_xy_res;
radius_boundingsphere_nucsize = radius_boundingsphere_nucsize/input_xy_res;


%everything below deals with raw values (timepoints, units in the representation)





featurevector = [];

for time_index = 1:length(times)
    
    disp(['time = ',num2str(time_index)])
    
    t = times(time_index);
    timepoint_featurevector = [];

    %aux

    %points from embryo--non-temporal features work within one time step
    points = embryo(t).finalpoints;

    %num points
    n = size(points,1);

    %normals, calculated with k_norm nns
    [~,normals,residuals] = cohereNormals(points,radius_boundingsphere_residual,3); % 
    %NOTE last argument of cohere normals is c = interconnectivity of graph


    for i=1:n
        
        disp(i)

        cell_featurevector = [];
        
        %skip if cell is a false
        if startsWith(embryo(t).cellnames{i},'Nuc')
            continue
        end
        

        %%
        %Cell Cycle

        %if any of the cell cycle metrics is requested, compute
        if sum(indicators(1:6)>0)

            %cell cycle metrics of current cell
            [beginning_of_current_cell_cycle,end_of_current_cell_cycle,duration_of_current_cell_cycle,time_until_next_division,time_since_last_division,...
                fraction_of_current_cycle_elapsed,time_of_next_division] = cellCycleMetrics(t,i,embryo);
            %leave out time of next division; redundant

            %add desired quantities to cell_featurevector
            if indicators(ordered_feats=="beginning_of_current_cell_cycle")
                cell_featurevector = [cell_featurevector,beginning_of_current_cell_cycle];
            end
            if indicators(ordered_feats=="end_of_current_cell_cycle")
                cell_featurevector = [cell_featurevector,end_of_current_cell_cycle];
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
                similarity_to_cousins_cell_cycle_cell_cycle = 1;
            else %assuming it's never zero
                similarity_to_cousins_cell_cycle_cell_cycle = cousin_metrics(end,3)/(mean(cousin_metrics(1:end-1,3))); %3 is location of cell cycle
                %ratio between cell i's cell cycle duration and cousins' cell
                %cycle durations
            end

            if size(ancestor_metrics,1) == 1
                similarity_to_ancestors_cell_cycle = 1;
            else %assuming it's never zero
                similarity_to_ancestors_cell_cycle = ancestor_metrics(end,3)/(mean(ancestor_metrics(1:end-1,3))); %3 is location of cell cycle
                %ratio between cell i's cell cycle duration and ancestors' cell
                %cycle durations
            end


            %add desired quantities to cell_featurevector
            if indicators(ordered_feats=="similarity_to_ancestors_cell_cycle")
                cell_featurevector = [cell_featurevector,similarity_to_ancestors_cell_cycle];
            end
            if indicators(ordered_feats=="similarity_to_cousins_cell_cycle")
                cell_featurevector = [cell_featurevector,similarity_to_cousins_cell_cycle_cell_cycle];
            end

        end

        %%
        %Cell density

        if sum(indicators(9:12)>0)

            %hyperparamters
            r = radius_boundingsphere_celldensity; 
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
        %Migration paths **RENAME vec_to_position_x_mins_ago...evc5 including the strings

        if sum(indicators(17:21)>0)

            %hyperparameter mins (each feature could have its own time window)
            tpts_back = timewindow_migrationpaths;%minutes
            
            %artifact from when parameter was input in minutes and
            %conversion was hardcoded - saving in case bug later
            %convert
            %mins = ceil(mins*305/470); %305 timesteps before twitching/470 mins (on average)


            %vec connecting cell to where it was x mins ago
            [past_time,past_ind,~,~] = backTraverse(t,i,'timesteps',tpts_back,embryo);
            vec_to_position_x_mins_ago = embryo(past_time).finalpoints(past_ind,:) - embryo(t).finalpoints(i,:);

            %vec connecting cell to where it will be x mins in future (or end of dataset)
            [future_pts] = forwardTraverse(t,[i],'timesteps',t+tpts_back,embryo);

            positions = [];
            
            
            
            for j = 1:length(future_pts)
                positions = [positions;embryo(min(t+tpts_back,length(embryo)-1)).finalpoints(future_pts(j),:)];
            end
            %in case cell/lineage dies
            if isempty(positions)
                vec_to_position_x_mins_in_future = [0,0,0];
            else
                positions = mean(positions,1);
                vec_to_position_x_mins_in_future = positions - embryo(t).finalpoints(i,:);
            end



            %vec connecting cell to site of last division
            [past_division_time,past_division_ind,~,~] = backTraverse(t,i,'generations',1,embryo);
            vec_to_site_of_last_division = embryo(past_division_time).finalpoints(past_division_ind,:) - embryo(t).finalpoints(i,:);

            %vec connecting cell to site of next division (0 if cell dies)
            timept = t;
            indpt = i;
            vec_to_site_of_next_division = [];
            if embryo(timept).suc(indpt,1) == -1
                vec_to_site_of_next_division = [0,0,0];
            else
                while embryo(timept).suc(indpt,2) == -1
                    indpt = embryo(timept).suc(indpt,1);
                    timept = timept + 1;
                    %if at end of embryo abort
                    if timept >= length(embryo)-1
                        vec_to_site_of_next_division = [0,0,0];
                        break
                    end
                    %abort if cell dies (at current or next time)
                    if embryo(timept).suc(indpt,1) == -1
                        vec_to_site_of_next_division = [0,0,0];
                        break
                    end
                end
            end

            if isempty(vec_to_site_of_next_division)
                vec_to_site_of_next_division = embryo(timept).finalpoints(indpt,:) - embryo(t).finalpoints(i,:);
            end


            %vec connecting site of last division to site of next division
            %uses previously defined quantities

            %handle having reached the end of embryo
            if timept > length(embryo)-1
                displacement_vec_over_cell_cycle = [0,0,0];
            else
                displacement_vec_over_cell_cycle = embryo(timept).finalpoints(indpt,:) - embryo(past_division_time).finalpoints(past_division_ind,:);
            end

            %add desired quantities to cell_featurevector
            if indicators(ordered_feats=="vec_to_position_x_mins_ago")
                cell_featurevector = [cell_featurevector,vec_to_position_x_mins_ago];
            end
            if indicators(ordered_feats=="vec_to_position_x_mins_in_future")
                cell_featurevector = [cell_featurevector,vec_to_position_x_mins_in_future];
            end
            if indicators(ordered_feats=="vec_to_site_of_last_division")
                cell_featurevector = [cell_featurevector,vec_to_site_of_last_division];
            end
            if indicators(ordered_feats=="vec_to_site_of_next_division")
                cell_featurevector = [cell_featurevector,vec_to_site_of_next_division];
            end
            if indicators(ordered_feats=="displacement_vec_over_cell_cycle")
                cell_featurevector = [cell_featurevector,displacement_vec_over_cell_cycle];
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
            timepoints_back = timewindow_migrationcongruency;
            
            %artifact from when parameter was input in minutes and
            %conversion was hardcoded - saving in case bug later
            %convert to timesteps
            %minutes_back = ceil(minutes_back*305/470);

            [previous_time,previous_ind,~,~] = backTraverse(t,i,'timesteps',timepoints_back,embryo);

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

            r_nuc = radius_boundingsphere_nucsize;

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

        if indicators(25)

            %developmental time (relative to 4 cell stage)***adjust to account for
            %temp/UV light--involves framerate (part of dataset)
            developmental_time = t;
            %accounting for 4cell stage time is done in normalization

            %add desired quantity to cell_featurevector (only feature for this category, so no need to check)
            cell_featurevector = [cell_featurevector,developmental_time];

        end


        %%

        %%% Feature Vector %%%

        timepoint_featurevector = [timepoint_featurevector;cell_featurevector];


    end
    
    featurevector = [featurevector;timepoint_featurevector];
    
    
end

%%% Normalize %%%

if normalize_bool

    %organize normalizations to be used
    kept_norms = ordered_norms(indicators==1);

    %make copy of featurevector
    normed_featurevector = featurevector;

    %aux
    
    %obtain last time at which there is an entry for timedata (movement starts)
    lasttime = length(embryo) - 1; %just a consequence of how the emb is loaded, assuming it is loaded up until movement starts

    %obtain max number of cells emb has before movement starts
    maxnumcells = size(embryo(lasttime).finalpoints,1);


    %compute normalizations
    for i = 1:size(indices_list,1)

        %get piece of normed_featurevector
        entries = normed_featurevector(:,indices_list(i,1):indices_list(i,2));

        %get norm type from hard-coded list corresponding to i'th entry of
        %feats_and_norms
        norm_type = kept_norms(kept_features == feats(i,1));

        if norm_type == "scale_4cell_to_end"
            %linearly scale all entries on [0,1] where 0 is when the emb is 4
            %cells and 1 is the end of the dataset (when the worm starts moving)
            entries = rescale(entries,0,1,'InputMin',fourcelltime,'InputMax',lasttime);

        elseif norm_type == "scale_0_to_max_num_cells"
            %linearly scale all entries on [0,1] where 0 is 4 cells (minimum number of cells if you begin counting time at the 4 cell stage) and
            %1 is the number of cells at end of the dataset (when the worm starts moving)
            entries = rescale(entries,0,1,'InputMin',4,'InputMax',maxnumcells);

        elseif norm_type == "column_zscore"
            %take zscore of column ( indices_list(i,1)==indices_list(i,2) )
            entries = zscore(entries);

        elseif norm_type == "all_zscore"
            %take zscore of all elements
            entries = zscore(entries,0,'all');

        elseif norm_type == "divide_avg_mag"
            %divide array elements by avg magitude of a vector in the array
            entries = avgMagnitude(entries);
        end

        %put newly normalized entries back in normed_featurevector
        normed_featurevector(:,indices_list(i,1):indices_list(i,2)) = entries;

    end
    
else    
    %if not computed
    normed_featurevector = [];
end


end
