%tuning hyperparameters mar 14

%for each case in cases, gives score = percent of 20 nns in feature space
%that are same type as a cell in question, averaged over all terminal cells

%requires preloading
%   embdat_stabilized
%   embeddingStruct
%   starttime
%   endtime
%   sampling
%   neighborparam


%%%


%inputs

%optimization jan25
%{
cases = load('optimization_cases_jan25.mat');
cases = cases.optimization_cases_jan25;
    %alternatively, the first column could be the features/normalizations and
    %the second column could be cell arrays containing the hyperparameter
    %values
%}


%feature evaluation mar 5
%{
%evaluate each feature separately

%start at 51 cells, tend = 305

%columns of output cell
%   1. feature computed
%   2. hyperparameter(s)
%   3. overall fate score
%   4. overall lineage score
%   5. percellfeature
%   6. fate_percellscores
%   7. lineage_percellscores

% NOTE: indices correspond to embeddingStruct loaded with :start at 51 cells, tend = 305
% NOTE: emb = '/Users/evanvietorisz/Documents/MATLAB/WT_embryos_Aug7/Decon_emb1_MGedits.zip';

% computed fractions using 20 nns
%}


%CURRENTLY REDUCES TO 15 DIMENSIONS MAY20


%cases = load('eval_cases_mar5_beforeresults.mat');
%cases = cases.eval_cases_mar5;

%cases = load('test_optimization_cases_mar14.mat');
%cases = cases.test_optimization_cases_mar14;

cases = load('optimization_cases_mar22.mat');
cases = cases.optimization_cases_mar22;

save_name = 'optimization_cases_may20_15dimensions_WITHRESULTS';



times = [starttime:sampling:endtime];
normalize_bool = true;
input_xy_res = .16;
    
number_neighbors = 20; %used to compute score (in feature space)


params_and_corresponding_feats = {"radius_boundingsphere_celldensity",["first_deriv_x";"first_deriv_y";"first_deriv_z";"laplacian"];
                                  "radius_boundingsphere_residual",   ["local_normal_vector";"residual"];
                                  "radius_boundingsphere_nucsize",    ["filtered_nuc_size"];
                                  "k_intercelldistance",              ["distances_to_k_nns";"signed_distances"];
                                  "degree_cousins",                   ["similarity_to_ancestors_cell_cycle";"similarity_to_cousins_cell_cycle"];
                                  "timewindow_migrationpaths",        ["vec_to_position_x_mins_ago";"vec_to_position_x_mins_in_future"];
                                  "timewindow_migrationcongruency",   ["proportion_of_present_nns_that_are_descendants_of_past_nns"];
                                  "k_migrationcongruency",            ["proportion_of_present_nns_that_are_descendants_of_past_nns"]};

params_and_corresponding_feats_wcommas = {"radius_boundingsphere_celldensity",["first_deriv_x","first_deriv_y","first_deriv_z","laplacian"];
                                  "radius_boundingsphere_residual",   ["local_normal_vector","residual"];
                                  "radius_boundingsphere_nucsize",    ["filtered_nuc_size"];
                                  "k_intercelldistance",              ["distances_to_k_nns","signed_distances"];
                                  "degree_cousins",                   ["similarity_to_ancestors_cell_cycle","similarity_to_cousins_cell_cycle"];
                                  "timewindow_migrationpaths",        ["vec_to_position_x_mins_ago","vec_to_position_x_mins_in_future"];
                                  "timewindow_migrationcongruency",   ["proportion_of_present_nns_that_are_descendants_of_past_nns"];
                                  "k_migrationcongruency",            ["proportion_of_present_nns_that_are_descendants_of_past_nns"]};






inds_to_keep = [1:length(embeddingStruct)]; %default to keeping all cells

%{
%Mar 8 creating dataset where you only take once of each cell per time

%%% PARED DOWN PERCELLFEATURE BELOW; UNDO THAT TO DO WORK ON PERCELLPERTIME
%%% FEATURESPACE

%obtain first indices of the each cell, eliminating repeated entries for
%the same cell at different time points
inds_to_keep = [];
cell_names = ["starter"];
for i = 1:length(embeddingStruct)
    if isempty(find(cell_names == convertCharsToStrings(embeddingStruct(i).name)))
        inds_to_keep = [inds_to_keep;i];
        cell_names = [cell_names;convertCharsToStrings(embeddingStruct(i).name)];
    end
    
end
%}



%backend

%misc
cases_output = cases;
embryo = embdat_stabilized;
     
%collect indices of each organ type
pharynxinds = [];
neuroninds = [];
hypoderminds = [];
bodywallmuscleinds = [];
gutinds = [];

for i = 1:length(embeddingStruct)
    if embeddingStruct(i).celltype == "pharynx"
        pharynxinds = [pharynxinds;i];
        
    elseif embeddingStruct(i).celltype == "neuron"
        neuroninds = [neuroninds;i];
        
    elseif embeddingStruct(i).celltype == "hypoderm"
        hypoderminds = [hypoderminds;i];
        
    elseif embeddingStruct(i).celltype == "bodywallmuscle"
        bodywallmuscleinds = [bodywallmuscleinds;i];
        
    elseif embeddingStruct(i).celltype == "gut"
        gutinds = [gutinds;i];

    end
end
    

%seed process by creating feats and param_values
feats = ["starter"];
param_values = cell(size(cases{1,2},1),2);
for i = 1:size(param_values,1)
    param_values{i,2} = 0;
end  % fill w zeros

%compute scores for all cases     
for i = 1:size(cases,1)
    
    %{ 
    %check if cases exists
    %see what is different about current cases from old cases
    %map parameters that changed to features affected (hardcoded cell)
    %get indices of the past features that changed and didn't change (feats
    %1,2,3 didn't change, but then feature 4 did change...etc)
    
    %*inside the try statement 
    %compute features that did change and are demanded in new feature set
    %assemble new featurevector
        %go down list of features demanded
        %if name array of features that didn't change from previous, put them
        %in (using indices previous feature computation)
        %if not, then go to name a array of features that did change from
        %previous (were computed) and put them in using indices output from
        %most recent feature computation
        %redefine variables to correspond to what is now the 'previous'
        %feature
    %} 

            
    new_feats = cases{i,1};
    new_param_values = cases{i,2};
    %new inds = blah?

    %map parameters that changed to features affected
    param_inds_that_changed = [new_param_values{:,2}] ~= [param_values{:,2}]; 
    feats_that_changed = [params_and_corresponding_feats_wcommas{param_inds_that_changed,2}]';

    %keep only feats that changed and were requested in new_feats
    %feats_that_changed = intersect(feats_that_changed,new_feats);


    %specify params
    %feats = cases{i,1};
    %param_values = cases{i,2};
    
    try
        
        %figure out what features need to be computed
        feats_to_compute = [];
        
        for j = 1:length(new_feats)%newly demanded features
            %if feature changed hyperparameter or was not present in
            %previous featuevector, compute it
            if ~isempty(find(feats_that_changed==new_feats(j))) || isempty(find(feats==new_feats(j)))
                feats_to_compute = [feats_to_compute;new_feats(j)];
            end
            
        end
        
                    
        %compute necessary features
        [computed_percellfeature,computed_normed_percellfeature,computed_indices_list] = computeFeatures4(embryo,times,feats_to_compute,new_param_values,normalize_bool,input_xy_res);
        
        
        
        %assemble featurevector
        new_percellfeature = [];
        new_normed_percellfeature = [];
        new_feature_sizes = [];
        
        for j = 1:length(new_feats)
            
            
            if ~isempty(find(feats_to_compute==new_feats(j))) %if feature was just computed, take from newly computed features
                
                beginning_ind = computed_indices_list(feats_to_compute==new_feats(j),1);
                end_ind = computed_indices_list(feats_to_compute==new_feats(j),2);
                
                new_percellfeature = [new_percellfeature,computed_percellfeature(:,beginning_ind:end_ind)];
                new_normed_percellfeature = [new_normed_percellfeature,computed_normed_percellfeature(:,beginning_ind:end_ind)];
                new_feature_sizes = [new_feature_sizes;(end_ind - beginning_ind + 1)];%repeated, take out if makes sense ^^^
                
            elseif ~isempty(find(feats==new_feats(j))) %otherwise, feature was present in previous featureset, so take form there
                
                beginning_ind = indices_list(feats==new_feats(j),1); %on first pass this statement won't get triggered, so indices_list isnt defined at outset
                end_ind = indices_list(feats==new_feats(j),2);
                
                new_percellfeature = [new_percellfeature,percellfeature(:,beginning_ind:end_ind)];
                new_normed_percellfeature = [new_normed_percellfeature,normed_percellfeature(:,beginning_ind:end_ind)];
                new_feature_sizes = [new_feature_sizes;(end_ind - beginning_ind + 1)];    
                
                
            else %otherwise there was an error
                disp('error')
            end

        end
        
        %redefine variables
        percellfeature = new_percellfeature;
        normed_percellfeature = new_normed_percellfeature;
        
        feats = new_feats;
        param_values = new_param_values;
        
        indices_list = zeros(size(new_feats,1),2);
        c = 1;
        for m = 1:size(new_feats,1)
            indices_list(m,1) = c;
            indices_list(m,2) = c + new_feature_sizes(m)-1;
            c = c + new_feature_sizes(m);
        end
        

        %**to look at per cell dataset only**
        adjusted_percellfeature = normed_percellfeature(inds_to_keep,:);
        adjusted_embeddingStruct = embeddingStruct(inds_to_keep); %make new embeddingStruct whose indices match pared down percellfeature
        
        cases_output{i,5} = adjusted_percellfeature;
        
        
    catch
        disp(['error in feature computation for case ',num2str(i)])
        cases_output{i,3} = "X";
        cases_output{i,4} = "X"; %indicates error occurred during feature computation
        cases_output{1,5} = [];
        continue
    end
    
    %*****lowering dimensions may 20*****
    
    %adjusted_percellfeature = run_umap(adjusted_percellfeature,'n_components',15,'n_neighbors',neighborparam);
    
    
    %compute scores

    lineage_percellscores = [];
    fate_percellscores = [];
    
    for j = 1:length(adjusted_percellfeature)
        
        %LINEAGE CONTINUITY

        %score is  percentage of nns that are same cell, progenitor, or
        %daughter
        
        %find neighbors
        [~,inds] = pdist2(adjusted_percellfeature,adjusted_percellfeature(j,:),'euclidean','Smallest',number_neighbors+1);
        inds = inds(2:end);

        %create list of names that qualify as same cell, progenitor, or
        %daughter
        
        %add same cell
        correct_names = [convertCharsToStrings(adjusted_embeddingStruct(j).name)];
        
        start_time = adjusted_embeddingStruct(j).captureFrame;
        start_ind = adjusted_embeddingStruct(j).captureIndex;
        
        %obtain name of progenitor
        
        [prog_time_pt,prog_ind_pt,~,~] = backTraverse(start_time,start_ind,'generations',1,embryo);
        prog_name = convertCharsToStrings(embryo(prog_time_pt).cellnames(prog_ind_pt));
        
        correct_names = [correct_names;prog_name];
        
        %obtain names of daughters
        [daughter_inds,daughter_time_pt] = forwardTraverse(start_time,start_ind,'generations',1,embryo);
        
        for l = 1:length(daughter_inds)
            correct_names = [correct_names;convertCharsToStrings(embryo(daughter_time_pt).cellnames{daughter_inds(l)})];
        end
        
        
        %tally how many neighbors are 'continuous'
        num_continuous = 0;
        for k = 1:number_neighbors

            if ~isempty(find(correct_names == adjusted_embeddingStruct(inds(k)).name))
                num_continuous = num_continuous+1;
            end

        end
        
        %find percent continuous, add to lineage_percellscores
        percent_continuous = num_continuous/number_neighbors;
        lineage_percellscores = [lineage_percellscores;percent_continuous];
        
        
        %TERMINAL FATE SIMILARITY
        if adjusted_embeddingStruct(j).celltype ~= "OTHER"
            %score is  percentage of nns that are same type as cell in
            %question
            
            %tally how many neighbors are of same type
            num_same_type = 0;
            for k = 1:number_neighbors
                
                if adjusted_embeddingStruct(inds(k)).celltype == adjusted_embeddingStruct(j).celltype
                    num_same_type = num_same_type+1;
                end
                
            end
            %percent of neighbors that are the same type
            percent_same_type = num_same_type/number_neighbors;
            
            %log score
            fate_percellscores = [fate_percellscores;percent_same_type];
            
        end
        
    end
    
    
    %log overall scores for that case
    
    fate_overall = mean(fate_percellscores);
    lineage_overall = mean(lineage_percellscores);
    
    cases_output{i,3} = fate_overall;
    cases_output{i,4} = lineage_overall;
    
    %log list of scores in case relevant later
    
    cases_output{i,6} = fate_percellscores;
    cases_output{i,7} = lineage_percellscores;
    
end

save(save_name,'cases_output','-v7.3')

%%
%making optimization cases mar 13

feats =      ["beginning_of_current_cell_cycle";
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

default_param_values = {"radius_boundingsphere_celldensity",20;
                        "radius_boundingsphere_residual",   20;
                        "radius_boundingsphere_nucsize",    20;
                        "k_intercelldistance",              50;
                        "degree_cousins",                   1;
                        "timewindow_migrationpaths",        .125;
                        "timewindow_migrationcongruency",   .125;
                        "k_migrationcongruency",            50};

param_ranges =         {"radius_boundingsphere_celldensity",[5:5:50];
                        "radius_boundingsphere_residual",   [5:5:50];
                        "radius_boundingsphere_nucsize",    [5:5:50];
                        "k_intercelldistance",              [15:5:50];
                        "degree_cousins",                   1;
                        "timewindow_migrationpaths",        [0.0250:0.0125:0.1125]; %10 mins to 45 mins in increments of 5 mins
                        "timewindow_migrationcongruency",   [0.0250:0.0125:0.1125];
                        "k_migrationcongruency",            [15:5:75]};


optimization_cases_mar13 = cell(1);

c = 1;

for i = 1:size(param_ranges,1)
    
    for j = 1:length(param_ranges{i,2})
        
        %alter the param values
        param_dummy = default_param_values;
        param_dummy{[param_dummy{:,1}]==param_ranges{i,1},2} = param_ranges{i,2}(j);
        
        
        optimization_cases_mar13{c,1} = feats;
        optimization_cases_mar13{c,2} = param_dummy;
        
        c = c+1;
        
    end
    
end



%%

%mar16 getting optimized embedding

default_param_values = {"radius_boundingsphere_celldensity",20;
                        "radius_boundingsphere_residual",   20;
                        "radius_boundingsphere_nucsize",    20;
                        "k_intercelldistance",              50;
                        "degree_cousins",                   1;
                        "timewindow_migrationpaths",        .125;
                        "timewindow_migrationcongruency",   .125;
                        "k_migrationcongruency",            50};
                    

                    
param_ranges =         {"radius_boundingsphere_celldensity",[5:5:50];
                        "radius_boundingsphere_residual",   [5:5:50];
                        "radius_boundingsphere_nucsize",    [5:5:50];
                        "k_intercelldistance",              [15:5:50];
                        "degree_cousins",                   1;
                        "timewindow_migrationpaths",        [0.0250:0.0125:0.1125]; %10 mins to 45 mins in increments of 5 mins
                        "timewindow_migrationcongruency",   [0.0250:0.0125:0.1125];
                        "k_migrationcongruency",            [15:5:75]};

num_trials_for_each_param = {"radius_boundingsphere_celldensity";
                        "radius_boundingsphere_residual";
                        "radius_boundingsphere_nucsize";
                        "k_intercelldistance";
                        "degree_cousins";
                        "timewindow_migrationpaths";
                        "timewindow_migrationcongruency";
                        "k_migrationcongruency"};
                    
for i = 1:size(param_ranges,1)
    num_trials_for_each_param{i,2} = length(param_ranges{i,2});
end


param_indices_list = zeros(size(num_trials_for_each_param,1),2);
c = 1;
for m = 1:size(num_trials_for_each_param,1)
    param_indices_list(m,1) = c;
    param_indices_list(m,2) = c + num_trials_for_each_param{m,2}-1;
    c = c + num_trials_for_each_param{m,2};
end
%%
%obtain best score for fate and lineage  

%make initial copies
fate_param_values = default_param_values;
lineage_param_values = default_param_values;

for i = 1:length(num_trials_for_each_param)
    
    cases_to_look_at = cases_output(param_indices_list(i,1):param_indices_list(i,2),:);
    
    [fate_best_score,fate_best_score_ind] = max([cases_to_look_at{:,3}]); %just happens to be 3rd column
    
    fate_best_param = cases_to_look_at{fate_best_score_ind,2}{i,2};     % grabs the value of the i'th parameter in the case 
                                                                        % that yielded the best score of the set (set is sweep 
                                                                        % over a given parameter)
    fate_param_values{i,2} = fate_best_param;
    
    [lineage_best_score,lineage_best_score_ind] = max([cases_to_look_at{:,4}]); %just happens to be 4th column
    lineage_best_param = cases_to_look_at{lineage_best_score_ind,2}{i,2};
    lineage_param_values{i,2} = lineage_best_param;
    
end

%%
% get plots (requires having things already loaded in workspace
% (neighborparam, celltype colors)

feats =      ["beginning_of_current_cell_cycle";
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
          
%make embeddings
[baseline_percellfeature_raw,baseline_percellfeature_normed,baseline_percellfeature_indices] = computeFeatures4(embryo,times,feats,default_param_values,normalize_bool,input_xy_res);
baseline_reduced = run_umap(baseline_percellfeature_normed,'n_neighbors',neighborparam);

[fate_percellfeature_raw,fate_percellfeature_normed,fate_percellfeature_indices] = computeFeatures4(embryo,times,feats,fate_param_values,normalize_bool,input_xy_res);
fate_reduced = run_umap(fate_percellfeature_normed,'n_neighbors',neighborparam);

[lineage_percellfeature_raw,lineage_percellfeature_normed,lineage_percellfeature_indices] = computeFeatures4(embryo,times,feats,lineage_param_values,normalize_bool,input_xy_res);
lineage_reduced = run_umap(lineage_percellfeature_normed,'n_neighbors',neighborparam);



%%
%save workspace
save('optimization_workspace_mar13.mat','-v7.3')


%%
%make plots (requires celltypecolors)

figure

subplot(1,3,1)
scatter(lineage_reduced(:,1),lineage_reduced(:,2),15,celltypecolors)
xlabel('UMAP-1')
ylabel('UMAP-2')
title('optimized for lineage continuity')
pbaspect([1,1,1])

subplot(1,3,2)
scatter(baseline_reduced(:,1),baseline_reduced(:,2),15,celltypecolors)
xlabel('UMAP-1')
ylabel('UMAP-2')
title('default parameters')
pbaspect([1,1,1])

subplot(1,3,3)
scatter(fate_reduced(:,1),fate_reduced(:,2),15,celltypecolors)
xlabel('UMAP-1')
ylabel('UMAP-2')
title('optimized for terminal fate clustering')
pbaspect([1,1,1])

sgtitle('after optimization, ALL features, mar13')

%%



for i = 1:length(num_trials_for_each_param)
    
    cases_to_look_at = optimization_cases_edited_mar13(param_indices_list(i,1):param_indices_list(i,2),:);
    
    
    
    optimization_cases_edited_mar13{1,2}{i,1}
    vals_swept = [param_ranges{i,2}]'
    fate_scores = round([cases_to_look_at{:,3}]',4)
    lineage_scores = round([cases_to_look_at{:,4}]',4)
    

end


%%

%%
%making optimization cases mar 22

feats =      ["beginning_of_current_cell_cycle";
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

default_param_values = {"radius_boundingsphere_celldensity",20;
                        "radius_boundingsphere_residual",   20;
                        "radius_boundingsphere_nucsize",    20;
                        "k_intercelldistance",              25;
                        "degree_cousins",                   1;
                        "timewindow_migrationpaths",        0.0750;
                        "timewindow_migrationcongruency",   0.0750;
                        "k_migrationcongruency",            10};

param_ranges =         {"radius_boundingsphere_celldensity",[5:5:50];
                        "radius_boundingsphere_residual",   [5:5:50];
                        "radius_boundingsphere_nucsize",    [5:5:50];
                        "k_intercelldistance",              [5:5:50];
                        "degree_cousins",                   1;
                        "timewindow_migrationpaths",        [0.0250:0.0125:0.1500]; %10 mins to 60 mins in increments of 5 mins
                        "timewindow_migrationcongruency",   [0.0250:0.0125:0.1500];
                        "k_migrationcongruency",            [5:5:25]};


optimization_cases_mar22 = cell(1);

c = 1;

for i = 1:size(param_ranges,1)
    
    for j = 1:length(param_ranges{i,2})
        
        %alter the param values
        param_dummy = default_param_values;
        param_dummy{[param_dummy{:,1}]==param_ranges{i,1},2} = param_ranges{i,2}(j);
        
        
        optimization_cases_mar22{c,1} = feats;
        optimization_cases_mar22{c,2} = param_dummy;
        
        c = c+1;
        
    end
    
end

%%
%mar 22 get optimized embedding
default_param_values = {"radius_boundingsphere_celldensity",20;
                        "radius_boundingsphere_residual",   20;
                        "radius_boundingsphere_nucsize",    20;
                        "k_intercelldistance",              25;
                        "degree_cousins",                   1;
                        "timewindow_migrationpaths",        0.0750;
                        "timewindow_migrationcongruency",   0.0750;
                        "k_migrationcongruency",            10};
                    
                    
param_ranges =         {"radius_boundingsphere_celldensity",[5:5:50];
                        "radius_boundingsphere_residual",   [5:5:50];
                        "radius_boundingsphere_nucsize",    [5:5:50];
                        "k_intercelldistance",              [5:5:50];
                        "degree_cousins",                   1;
                        "timewindow_migrationpaths",        [0.0250:0.0125:0.1500]; %10 mins to 60 mins in increments of 5 mins
                        "timewindow_migrationcongruency",   [0.0250:0.0125:0.1500];
                        "k_migrationcongruency",            [5:5:50]};

num_trials_for_each_param = param_ranges(:,1);
for i = 1:size(param_ranges,1)
    num_trials_for_each_param{i,2} = length(param_ranges{i,2});
end




param_indices_list = zeros(size(num_trials_for_each_param,1),2);
c = 1;
for m = 1:size(num_trials_for_each_param,1)
    param_indices_list(m,1) = c;
    param_indices_list(m,2) = c + num_trials_for_each_param{m,2}-1;
    c = c + num_trials_for_each_param{m,2};
end

%%
%obtain best score for fate and lineage  

%make initial copies
fate_param_values = default_param_values;
lineage_param_values = default_param_values;

for i = 1:length(num_trials_for_each_param)
    
    cases_to_look_at = cases_output(param_indices_list(i,1):param_indices_list(i,2),:);
    
    [fate_best_score,fate_best_score_ind] = max([cases_to_look_at{:,3}]); %just happens to be 3rd column
    
    fate_best_param = cases_to_look_at{fate_best_score_ind,2}{i,2};     % grabs the value of the i'th parameter in the case 
                                                                        % that yielded the best score of the set (set is sweep 
                                                                        % over a given parameter)
    fate_param_values{i,2} = fate_best_param;
    
    [lineage_best_score,lineage_best_score_ind] = max([cases_to_look_at{:,4}]); %just happens to be 4th column
    lineage_best_param = cases_to_look_at{lineage_best_score_ind,2}{i,2};
    lineage_param_values{i,2} = lineage_best_param;
    
end

%%
%mar 22
% get plots (requires having things already loaded in workspace
% (neighborparam, celltype colors)

feats =      ["beginning_of_current_cell_cycle";
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
          
%make embeddings
[baseline_percellfeature_raw,baseline_percellfeature_normed,baseline_percellfeature_indices] = computeFeatures4(embryo,times,feats,default_param_values,normalize_bool,input_xy_res);
baseline_reduced = run_umap(baseline_percellfeature_normed,'n_neighbors',neighborparam);

[fate_percellfeature_raw,fate_percellfeature_normed,fate_percellfeature_indices] = computeFeatures4(embryo,times,feats,fate_param_values,normalize_bool,input_xy_res);
fate_reduced = run_umap(fate_percellfeature_normed,'n_neighbors',neighborparam);

[lineage_percellfeature_raw,lineage_percellfeature_normed,lineage_percellfeature_indices] = computeFeatures4(embryo,times,feats,lineage_param_values,normalize_bool,input_xy_res);
lineage_reduced = run_umap(lineage_percellfeature_normed,'n_neighbors',neighborparam);

%%
%make plots (requires celltypecolors)

figure

subplot(1,3,1)
scatter(lineage_reduced(:,1),lineage_reduced(:,2),15,celltypecolors)
xlabel('UMAP-1')
ylabel('UMAP-2')
title('optimized for lineage continuity')
pbaspect([1,1,1])

subplot(1,3,2)
scatter(baseline_reduced(:,1),baseline_reduced(:,2),15,celltypecolors)
xlabel('UMAP-1')
ylabel('UMAP-2')
title('default parameters')
pbaspect([1,1,1])

subplot(1,3,3)
scatter(fate_reduced(:,1),fate_reduced(:,2),15,celltypecolors)
xlabel('UMAP-1')
ylabel('UMAP-2')
title('optimized for terminal fate clustering')
pbaspect([1,1,1])

sgtitle('after optimization, ALL features, mar22')

%opening optimization_WT_mar22_workspace brings you to here

%%
for i = 1:length(num_trials_for_each_param)
    
    cases_to_look_at = cases_output(param_indices_list(i,1):param_indices_list(i,2),:);
    
    
    
    optimization_cases_mar22{1,2}{i,1}
    vals_swept = [param_ranges{i,2}]'
    fate_scores = round([cases_to_look_at{:,3}]',4)
    lineage_scores = round([cases_to_look_at{:,4}]',4)
    
    

end

%%

%looking at lineages in WT

%obtain indices indices corresponding to each lineage

ABa_inds = [];
ABp_inds = [];
C_inds = [];
D_inds = [];
E_inds = [];
MS_inds = [];

for i = 1:length(embeddingStruct)
    if embeddingStruct(i).lineage == "ABa"
        ABa_inds = [ABa_inds;i];
    elseif embeddingStruct(i).lineage == "ABp"
        ABp_inds = [ABp_inds;i];
    elseif embeddingStruct(i).lineage == "C"
        C_inds = [C_inds;i];
    elseif embeddingStruct(i).lineage == "D"
        D_inds = [D_inds;i];
    elseif embeddingStruct(i).lineage == "E"
        E_inds = [E_inds;i];
    elseif embeddingStruct(i).lineage == "MS"
        MS_inds = [MS_inds;i];
    end
end

%looking at lineage results bc why not

pcf = lineage_percellfeature_normed;

ABa_pcf = pcf(ABa_inds,:);
ABp_pcf = pcf(ABp_inds,:);
C_pcf = pcf(C_inds,:);
D_pcf = pcf(D_inds,:);
E_pcf = pcf(E_inds,:);
MS_pcf = pcf(MS_inds,:);

ABa_reduced = run_umap(ABa_pcf,'n_neighbors',neighborparam);
ABp_reduced = run_umap(ABp_pcf,'n_neighbors',neighborparam);
C_reduced = run_umap(C_pcf,'n_neighbors',neighborparam);
D_reduced = run_umap(D_pcf,'n_neighbors',neighborparam);
E_reduced = run_umap(E_pcf,'n_neighbors',neighborparam);
MS_reduced = run_umap(MS_pcf,'n_neighbors',neighborparam);

ABa_celltypecolors = celltypecolors(ABa_inds,:);
ABp_celltypecolors = celltypecolors(ABp_inds,:);
C_celltypecolors = celltypecolors(C_inds,:);
D_celltypecolors = celltypecolors(D_inds,:);
E_celltypecolors = celltypecolors(E_inds,:);
MS_celltypecolors = celltypecolors(MS_inds,:);




figure

subplot(2,3,1)
scatter(ABa_reduced(:,1),ABa_reduced(:,2),20,ABa_celltypecolors)
title('ABa')
xlabel('UMAP-1')
ylabel('UMAP-2')

subplot(2,3,2)
scatter(ABp_reduced(:,1),ABp_reduced(:,2),20,ABp_celltypecolors)
title('ABp')
xlabel('UMAP-1')
ylabel('UMAP-2')

subplot(2,3,3)
scatter(C_reduced(:,1),C_reduced(:,2),20,C_celltypecolors)
title('C')
xlabel('UMAP-1')
ylabel('UMAP-2')

subplot(2,3,4)
scatter(D_reduced(:,1),D_reduced(:,2),20,D_celltypecolors)
title('D')
xlabel('UMAP-1')
ylabel('UMAP-2')

subplot(2,3,5)
scatter(E_reduced(:,1),E_reduced(:,2),20,E_celltypecolors)
title('E')
xlabel('UMAP-1')
ylabel('UMAP-2')

subplot(2,3,6)
scatter(MS_reduced(:,1),MS_reduced(:,2),20,MS_celltypecolors)
title('MS')
xlabel('UMAP-1')
ylabel('UMAP-2')

tt = 'lineages embedded separately | using lineage-optimized feature space';
tt = [tt newline 'Pharynx (g), Neuron (b), Body Wall Muscle (r), Hypoderm (f), Gut (c)']; 
sgtitle(tt)





%%

%getting time colors


ABa_timecolors = timecolors(ABa_inds,:);
ABp_timecolors = timecolors(ABp_inds,:);
C_timecolors = timecolors(C_inds,:);
D_timecolors = timecolors(D_inds,:);
E_timecolors = timecolors(E_inds,:);
MS_timecolors = timecolors(MS_inds,:);


figure

subplot(2,3,1)
scatter(ABa_reduced(:,1),ABa_reduced(:,2),20,ABa_timecolors)
title('ABa')
xlabel('UMAP-1')
ylabel('UMAP-2')

subplot(2,3,2)
scatter(ABp_reduced(:,1),ABp_reduced(:,2),20,ABp_timecolors)
title('ABp')
xlabel('UMAP-1')
ylabel('UMAP-2')

subplot(2,3,3)
scatter(C_reduced(:,1),C_reduced(:,2),20,C_timecolors)
title('C')
xlabel('UMAP-1')
ylabel('UMAP-2')

subplot(2,3,4)
scatter(D_reduced(:,1),D_reduced(:,2),20,D_timecolors)
title('D')
xlabel('UMAP-1')
ylabel('UMAP-2')

subplot(2,3,5)
scatter(E_reduced(:,1),E_reduced(:,2),20,E_timecolors)
title('E')
xlabel('UMAP-1')
ylabel('UMAP-2')

subplot(2,3,6)
scatter(MS_reduced(:,1),MS_reduced(:,2),20,MS_timecolors)
title('MS')
xlabel('UMAP-1')
ylabel('UMAP-2')

tt = 'lineages embedded separately | using lineage-optimized feature space';
tt = [tt newline 'time: green -> red (51 cells to 571 cells)']; 
sgtitle(tt)


%%
%scoring the WT

[BASE_fate_overall,BASE_lineage_overall,~,~] = fateAndLineageScores(embryo,baseline_percellfeature_normed,embeddingStruct,20);
[FATE_fate_overall,FATE_lineage_overall,~,~] = fateAndLineageScores(embryo,fate_percellfeature_normed,embeddingStruct,20);
[LINEAGE_fate_overall,LINEAGE_lineage_overall,~,~] = fateAndLineageScores(embryo,lineage_percellfeature_normed,embeddingStruct,20);

%%

%lineage comparisons
%make a plot where the lineages are embedded together in differnt point
%styles, colored as they would be by themselves

ABa_inds = [];
ABp_inds = [];
C_inds = [];
D_inds = [];
E_inds = [];
MS_inds = [];

for i = 1:length(embeddingStruct)
    if embeddingStruct(i).lineage == "ABa"
        ABa_inds = [ABa_inds;i];
    elseif embeddingStruct(i).lineage == "ABp"
        ABp_inds = [ABp_inds;i];
    elseif embeddingStruct(i).lineage == "C"
        C_inds = [C_inds;i];
    elseif embeddingStruct(i).lineage == "D"
        D_inds = [D_inds;i];
    elseif embeddingStruct(i).lineage == "E"
        E_inds = [E_inds;i];
    elseif embeddingStruct(i).lineage == "MS"
        MS_inds = [MS_inds;i];
    end
end

%looking at lineage results bc why not
pcf = lineage_percellfeature_normed;



%define two lineages to put in
lin1_pcf  = pcf(C_inds,:);
lin1_colors = celltypecolors(C_inds,:);
lin1_colors_time = timecolors(C_inds,:);

lin2_pcf =  pcf(MS_inds,:);
lin2_colors = celltypecolors(MS_inds,:);
lin2_colors_time = timecolors(MS_inds,:);

lin1_name = 'C';
lin2_name = 'MS';




reduced_feat = run_umap([lin1_pcf;lin2_pcf],'n_neighbors',neighborparam);
reduced_feat_lin1 = reduced_feat(1:size(lin1_pcf,1),:);
reduced_feat_lin2 = reduced_feat(size(lin1_pcf,1)+1:end,:);



figure

subplot(1,2,1)
hold on
scatter(reduced_feat_lin1(:,1),reduced_feat_lin1(:,2),20,lin1_colors)
scatter(reduced_feat_lin2(:,1),reduced_feat_lin2(:,2),20,lin2_colors,'filled','d','MarkerEdgeColor',[0,0,0])
hold off
xlabel('UMAP-1')
ylabel('UMAP-2')
title('Pharynx (g), Neuron (b), Body Wall Muscle (r), Hypoderm (f), Gut (c)')

subplot(1,2,2)
hold on
scatter(reduced_feat_lin1(:,1),reduced_feat_lin1(:,2),20,lin1_colors_time)
scatter(reduced_feat_lin2(:,1),reduced_feat_lin2(:,2),20,lin2_colors_time,'filled','d','MarkerEdgeColor',[0,0,0])
hold off
xlabel('UMAP-1')
ylabel('UMAP-2')
title('time: green -> red (51 cells to 571 cells)')

sgtitle(['Co-decomposition: ',lin1_name,' (circles), ',lin2_name,' (diamonds)'])



%%

%lineages from the whole embedding by separated out

pcf = lineage_percellfeature_normed;

reduced_feat = run_umap(pcf,'n_neighbors',neighborparam);

%lineage comparisons
%make a plot where the lineages are embedded together in differnt point
%styles, colored as they would be by themselves

ABa_inds = [];
ABp_inds = [];
C_inds = [];
D_inds = [];
E_inds = [];
MS_inds = [];

for i = 1:length(embeddingStruct)
    if embeddingStruct(i).lineage == "ABa"
        ABa_inds = [ABa_inds;i];
    elseif embeddingStruct(i).lineage == "ABp"
        ABp_inds = [ABp_inds;i];
    elseif embeddingStruct(i).lineage == "C"
        C_inds = [C_inds;i];
    elseif embeddingStruct(i).lineage == "D"
        D_inds = [D_inds;i];
    elseif embeddingStruct(i).lineage == "E"
        E_inds = [E_inds;i];
    elseif embeddingStruct(i).lineage == "MS"
        MS_inds = [MS_inds;i];
    end
end


%%

%define two lineages to put in
lin1_reduced  = reduced_feat(MS_inds,:);
lin1_colors = celltypecolors(MS_inds,:);
lin1_colors_time = timecolors(MS_inds,:);
reduced_feat_lin1 = reduced_feat(MS_inds,:);

lin2_pcf =  reduced_feat(D_inds,:);
lin2_colors = celltypecolors(D_inds,:);
lin2_colors_time = timecolors(D_inds,:);
reduced_feat_lin2 = reduced_feat(D_inds,:);

lin1_name = 'MS';
lin2_name = 'D';



figure

subplot(1,2,1)
hold on
scatter(reduced_feat_lin1(:,1),reduced_feat_lin1(:,2),20,lin1_colors)
scatter(reduced_feat_lin2(:,1),reduced_feat_lin2(:,2),6,lin2_colors,'filled','MarkerEdgeColor',[0,0,0])
hold off
xlabel('UMAP-1')
ylabel('UMAP-2')
title('Pharynx (g), Neuron (b), Body Wall Muscle (r), Hypoderm (f), Gut (c)')

subplot(1,2,2)
hold on
scatter(reduced_feat_lin1(:,1),reduced_feat_lin1(:,2),20,lin1_colors_time)
scatter(reduced_feat_lin2(:,1),reduced_feat_lin2(:,2),6,lin2_colors_time,'filled','MarkerEdgeColor',[0,0,0])
hold off
xlabel('UMAP-1')
ylabel('UMAP-2')
title('time: green -> red (51 cells to 571 cells)')

sgtitle(['lineages taken from embryo embedding: ',lin1_name,' (rings), ',lin2_name,' (dots)'])


%%

%figure to show all lineages separately

%get bounds

xbounds = [min(reduced_feat(:,1)-.1),max(reduced_feat(:,1))+.1];
ybounds = [min(reduced_feat(:,2)-.1),max(reduced_feat(:,2))+.1];


figure
subplot(2,4,1)
scatter(reduced_feat(ABa_inds,1),reduced_feat(ABa_inds,2),20,celltypecolors(ABa_inds,:))
xlabel('UMAP-1')
ylabel('UMAP-2')
title('ABa')
xlim(xbounds)
ylim(ybounds)

subplot(2,4,2)
scatter(reduced_feat(ABp_inds,1),reduced_feat(ABp_inds,2),20,celltypecolors(ABp_inds,:))
xlabel('UMAP-1')
ylabel('UMAP-2')
title('ABp')
xlim(xbounds)
ylim(ybounds)

subplot(2,4,3)
scatter(reduced_feat(C_inds,1),reduced_feat(C_inds,2),20,celltypecolors(C_inds,:))
xlabel('UMAP-1')
ylabel('UMAP-2')
title('C')
xlim(xbounds)
ylim(ybounds)

subplot(2,4,5)
scatter(reduced_feat(D_inds,1),reduced_feat(D_inds,2),20,celltypecolors(D_inds,:))
xlabel('UMAP-1')
ylabel('UMAP-2')
title('D')
xlim(xbounds)
ylim(ybounds)

subplot(2,4,6)
scatter(reduced_feat(E_inds,1),reduced_feat(E_inds,2),20,celltypecolors(E_inds,:))
xlabel('UMAP-1')
ylabel('UMAP-2')
title('E')
xlim(xbounds)
ylim(ybounds)

subplot(2,4,7)
scatter(reduced_feat(MS_inds,1),reduced_feat(MS_inds,2),20,celltypecolors(MS_inds,:))
xlabel('UMAP-1')
ylabel('UMAP-2')
title('MS')
xlim(xbounds)
ylim(ybounds)

%two full embeddings

subplot(2,4,4)
scatter(reduced_feat(:,1),reduced_feat(:,2),20,celltypecolors)
xlabel('UMAP-1')
ylabel('UMAP-2')
title('All Lineages')
xlim(xbounds)
ylim(ybounds)

subplot(2,4,8)
scatter(reduced_feat(:,1),reduced_feat(:,2),20,timecolors)
xlabel('UMAP-1')
ylabel('UMAP-2')
xlim(xbounds)
ylim(ybounds)
title('time: green -> red (51 cells to 571 cells)')

tt = 'individual lineages taken from entire embedding';
tt = [tt newline 'Pharynx (g), Neuron (b), Body Wall Muscle (r), Hypoderm (f), Gut (c)'];
sgtitle(tt)
%%

%tracing individual lineages


%good ABs
%35894

%good MSs
%35903


%define two cells

cell1_structIndex = 35894;
cell1_color = [1,0,0];

cell2_structIndex = 35903;
cell2_color = [0,1,0];

overlap_color = [.5,0,1];

%color colors approproately
custom_colors = repmat([.5,.5,.5],length(embeddingStruct),1);


disp(embeddingStruct(cell1_structIndex).name)
disp(embeddingStruct(cell2_structIndex).name)
disp(' ')



for i = length(embeddingStruct):-1:1
    
    if startsWith(embeddingStruct(cell1_structIndex).name,embeddingStruct(i).name)
        custom_colors(i,:) = cell1_color;
        
        
        
        
        if startsWith(embeddingStruct(cell2_structIndex).name,embeddingStruct(i).name)
            custom_colors(i,:) = overlap_color;
        end
 
    elseif startsWith(embeddingStruct(cell2_structIndex).name,embeddingStruct(i).name)
         custom_colors(i,:) = cell2_color;
         
         
    end
    
end


figure
hold on
scatter(reduced_feat(ABa_inds,1),reduced_feat(ABa_inds,2),20,custom_colors(ABa_inds,:))

scatter(reduced_feat(MS_inds,1),reduced_feat(MS_inds,2),10,custom_colors(MS_inds,:),'filled','MarkerEdgeColor',[0,0,0])
hold off


%%


%coords in embedding %%%currently mutant

reduced_to_look_at = reduced_all(WT_indices,:);
struct_to_look_at = embeddingStruct_truncated;

cols = repmat([.25,.25,.25],size(reduced_to_look_at,1),1);


xmin = .74;
xmax =4;

ymin =2.5;
ymax = 4;

for i = 1:size(reduced_to_look_at,1)
    
    if reduced_to_look_at(i,1) >= xmin && reduced_to_look_at(i,1) <= xmax &&  reduced_to_look_at(i,2) >= ymin && reduced_to_look_at(i,2) <= ymax
        disp(struct_to_look_at(i).name)
        cols(i,:) = [1,0,0];
    end
    
    
end
%%
scatter(reduced_to_look_at(:,1),reduced_to_look_at(:,2),20,cols)
















