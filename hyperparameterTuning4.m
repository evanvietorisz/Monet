%tuning hyperparameters Jun 29
%   uses scoring metrics where lineage is constrained in time and fate is
%   done on per cell dataset (one instance of each cellID instead of
%   multiple)

%different from hyper3 because the best scoring parameter is updated in
%all remaining trials


%{

general layout

makes a cell called optimization cases that contains differnet parameter
combinations, saves it

then, loads those parameter combinations, computes the feature space using
them, scores it, and saves the scores in a new, similarly shaped cell
called optimizationcases_WITHRESULTS.

one can then open that data to look at which parameter combinations
performed the best and recover the feature space that went along with them

optimizes by setting all params to default values, then optimizes each one
in sequence over the range specified in param_ranges. after each sweep,
picks the best value and sets it to that for sweeps of subsequent params




%}


%%

%%% generate optimization cases %%%

%aux--note: additional useful scripts in the 'hyperparameterTuning2' underneath
%'make featurevector mar 22' section

%make cases jun 10

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

          %starting values for each param
default_param_values = {"radius_boundingsphere_celldensity",15;
                        "radius_boundingsphere_residual",   15;
                        "radius_boundingsphere_nucsize",    15;
                        "k_intercelldistance",              25;
                        "degree_cousins",                   1;
                        "timewindow_migrationpaths",        0.0750;
                        "timewindow_migrationcongruency",   0.0750;
                        "k_migrationcongruency",            15};

param_ranges =         {"radius_boundingsphere_celldensity",[5:5:25];
                        "radius_boundingsphere_residual",   [5:5:25];
                        "radius_boundingsphere_nucsize",    [5:5:25];
                        "k_intercelldistance",              [5:5:50];
                        "degree_cousins",                   1;
                        "timewindow_migrationpaths",        [0.0250:0.0125:0.1500]; %10 mins to 60 mins in increments of 5 mins
                        "timewindow_migrationcongruency",   [0.0250:0.0125:0.1500];
                        "k_migrationcongruency",            [5:5:25]};

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
                    
optimization_cases_jun29 = cell(1);

c = 1;

for i = 1:size(param_ranges,1)
    
    for j = 1:length(param_ranges{i,2})
        
        %alter the param values
        param_dummy = default_param_values;
        param_dummy{[param_dummy{:,1}]==param_ranges{i,1},2} = param_ranges{i,2}(j);
        
        
        optimization_cases_jun29{c,1} = feats;
        optimization_cases_jun29{c,2} = param_dummy;
        
        c = c+1;
        
    end
    
end

save('optimization_cases_jun29.mat')


%%


%%% uploading files, preprocessing %%%

%copied from celeganstest jun 10

%load embryo
emb = '/Users/evanvietorisz/Documents/MATLAB/WT_embryos_Aug7/Decon_emb1_MGedits.zip';
input_xy_res = .16; % resolution of xy axes in microns
anisotropy = 1;

tstart=1;
tend=305; %cutoff should be when worm starts moving around (determined manually)

centering_time = tend; %last trustable timestep in dataset--should be same as tend
sampling=2;
endtime=tend-1; %should be one less than last timepoint loaded because of implementation

starting_num_cells = 51; %number of cells at which to start computing features--must be one greater than k

normalize_bool = true;

neighborparam = 30; %perplexity in tsne and n_neighbors in umap

time_window = 7.5; %number of minutes forward and back to look when evaluating lineage scores



%load and prepare data
templocation='temp_unzip/';
%unzip zipfile to temp file
if ~exist(emb,'file')
    errors.zipfilemissing=1;
    return
end
disp('the file exists')
try
    unzip(emb,templocation);
catch exception
    errors.zipfilecorrupted=1;
    return
end
disp('the file was unzipped')

%load embryo struct and cells struct
[cells,embdat] = loadcells_unnamed(templocation,tend,anisotropy,false);
rmdir(templocation,'s');


%align input to reference (reference is NOT internally aligned)
%note the .16 is  scale of the reference and is just a property of the data
load('ref_emb.mat')
embdat_stabilized = coalignNamedEmbryosPerTime2(ref_emb,1,size(ref_emb,2),.16,embdat,1,tend,anisotropy,input_xy_res,centering_time);

%note: anisotropy appears not to be necessary for the alignment-- no diff
%difference
disp('the input was transformed')

%update cells object
cells_stabilized = parseCellsFromEmb(embdat_stabilized,tend);

%backup
cells_backup=cells;
cells=cells_stabilized;
disp('the data was prepared')



%start at specified number of cells***abstract this away
starttime = 1;
while size(embdat_stabilized(starttime).finalpoints,1) < starting_num_cells
    starttime = starttime + 1;
end

%embeddingStruct
celltypes = ["pharynx","neuron","bodywallmuscle","hypoderm","gut"];
lineages = ["ABa","ABp","C","D","E","MS","P"];

times = [starttime:sampling:endtime];

embeddingStruct = embeddingStructMaker(embdat_stabilized,starttime,sampling,endtime,lineages);

%remove any debris/false nucs
nuc_inds = [];
for i = 1:length(embeddingStruct)
    if startsWith(embeddingStruct(i).name,'Nuc')
        nuc_inds = [nuc_inds;i];
    end
end
embeddingStruct(nuc_inds) = [];

numPointsInEmbedding = size(embeddingStruct,2);


%%

%%% optimization %%%

%inputs

cases = load('optimization_cases_jun29.mat');
cases = cases.optimization_cases_jun29;

save_name = 'optimization_cases_jun29_WITHRESULTS';
    
number_neighbors = 20; %used to compute score (in feature space)

params_and_corresponding_feats_wcommas = {"radius_boundingsphere_celldensity",["first_deriv_x","first_deriv_y","first_deriv_z","laplacian"];
                                  "radius_boundingsphere_residual",   ["local_normal_vector","residual"];
                                  "radius_boundingsphere_nucsize",    ["filtered_nuc_size"];
                                  "k_intercelldistance",              ["distances_to_k_nns","signed_distances"];
                                  "degree_cousins",                   ["similarity_to_ancestors_cell_cycle","similarity_to_cousins_cell_cycle"];
                                  "timewindow_migrationpaths",        ["vec_to_position_x_mins_ago","vec_to_position_x_mins_in_future"];
                                  "timewindow_migrationcongruency",   ["proportion_of_present_nns_that_are_descendants_of_past_nns"];
                                  "k_migrationcongruency",            ["proportion_of_present_nns_that_are_descendants_of_past_nns"]};



inds_to_keep = [1:length(embeddingStruct)]; %default to keeping all cells

%backend

%misc
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
end 


%%%changes start here jun 29

cases_output = cell(size(cases,1),6);

optimize_for = 'lineage';

%compute param indices 
param_indices_list = zeros(size(num_trials_for_each_param,1),2);
c = 1;
for m = 1:size(num_trials_for_each_param,1)
    param_indices_list(m,1) = c;
    param_indices_list(m,2) = c + num_trials_for_each_param{m,2}-1;
    c = c + num_trials_for_each_param{m,2};
end

%compute all cases, update params as you go
num_of_params = size(default_param_values,1); %for each parameter
for pdex = 1:num_of_params
    
    param_cases = cases(param_indices_list(pdex,1):param_indices_list(pdex,2),:);%the piece of cases taken out from
    
    %do the trials, compute scores
    for i = 1:size(param_cases,1)

        new_feats = param_cases{i,1};
        new_param_values = param_cases{i,2};

        %map parameters that changed to features affected
        param_inds_that_changed = [new_param_values{:,2}] ~= [param_values{:,2}]; 
        feats_that_changed = [params_and_corresponding_feats_wcommas{param_inds_that_changed,2}]';


        %compute new featurevector

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


        %ability to select a subset of cells in the emb if you want
        adjusted_percellfeature = normed_percellfeature(inds_to_keep,:);
        adjusted_embeddingStruct = embeddingStruct(inds_to_keep); %make new embeddingStruct whose indices match pared down percellfeature

        param_cases{i,5} = adjusted_percellfeature;


        %lowering dimensions for consistency on different featurevector sizes
        lowD_percellfeature = run_umap(adjusted_percellfeature,'n_components',15,'n_neighbors',neighborparam);

        param_cases{i,6} = lowD_percellfeature;

        %compute scores

        %make condensed version for fate scoring
        [condensed_lowD_percellfeature,condensed_embeddingStruct] = condenseToOneEntryPerCell(lowD_percellfeature,adjusted_embeddingStruct,'median_timepoint_entry');

        %score, log
        [fate_overall,lineage_overall]  = computeScore3(lowD_percellfeature,adjusted_embeddingStruct,embryo,number_neighbors,time_window,condensed_lowD_percellfeature,condensed_embeddingStruct,starttime,sampling,endtime);

        %log overall scores for that case

        param_cases{i,3} = fate_overall;
        param_cases{i,4} = lineage_overall;

    end
    
    %find the best trial from that sweep
    switch optimize_for
        case 'lineage'
            struct_index_to_get_score_from = 4;
        case 'fate'
            struct_index_to_get_score_from = 4;
    end
    [best_score_val, best_score_ind] = max([param_cases{:,4}]);%4 makes it lineage score, 3 is fate
    best_param_value = param_cases{best_score_ind,2}{pdex,2};
    
    %put the new trials from the sweep back in the cases cell
    cases_output(param_indices_list(pdex,1):param_indices_list(pdex,2),:) = param_cases;
    
    %put the optimized value in the remaining cases
    if pdex ~= num_of_params
        for m = param_indices_list(pdex,2)+1:size(cases,1)%all cases in subsequent sweeps
            cases{m,2}{pdex,2} = best_param_value;
        end
    end
    
    
end

%obtain the best param combination by looking at the final sweep
best_params_case_entry = param_cases(1,:);

[best_score_val, best_score_ind] = max([param_cases{:,4}]);
best_param_value = param_cases{best_score_ind,2}{pdex,2};

best_params_case_entry{2}{pdex,2} = best_param_value;


save(save_name,'cases_output','-v7.3')%the best_params_case_entry saved with this is not correct!!!!only the param values are correct


%%
%compute, score best case entry june 30

optimized_param_values = best_params_case_entry{2};

%make embeddings
[~,baseline_percellfeature_normed,~] = computeFeatures4(embryo,times,feats,default_param_values,normalize_bool,input_xy_res);
baseline_reduced = run_umap(baseline_percellfeature_normed,'n_neighbors',neighborparam);

[~,optimized_percellfeature_normed,~] = computeFeatures4(embryo,times,feats,optimized_param_values,normalize_bool,input_xy_res);
optimized_reduced = run_umap(optimized_percellfeature_normed,'n_neighbors',neighborparam);

save('optimization_jun29_results_files.mat','baseline_percellfeature_normed','optimized_percellfeature_normed','baseline_reduced','optimized_reduced','default_param_values','optimized_param_values','-v7.3')

%score
baseline_lowD_percellfeature = run_umap(baseline_percellfeature_normed,'n_components',15,'n_neighbors',neighborparam);
[baseline_condensed_lowD_percellfeature,baseline_condensed_embeddingStruct] = condenseToOneEntryPerCell(baseline_lowD_percellfeature,adjusted_embeddingStruct,'median_timepoint_entry');
[baseline_fate_score, baseline_lineage_score] = computeScore3(baseline_lowD_percellfeature,adjusted_embeddingStruct,embryo,number_neighbors,time_window,baseline_condensed_lowD_percellfeature,baseline_condensed_embeddingStruct,starttime,sampling,endtime);

optimized_lowD_percellfeature = run_umap(optimized_percellfeature_normed,'n_components',15,'n_neighbors',neighborparam);
[optimized_condensed_lowD_percellfeature,optimized_condensed_embeddingStruct] = condenseToOneEntryPerCell(optimized_lowD_percellfeature,adjusted_embeddingStruct,'median_timepoint_entry');

[optimized_fate_score, optimized_lineage_score] = computeScore3(optimized_lowD_percellfeature,adjusted_embeddingStruct,embryo,number_neighbors,time_window,optimized_condensed_lowD_percellfeature,optimized_condensed_embeddingStruct,starttime,sampling,endtime);

save('optimized_scores_jun29.mat','baseline_fate_score','baseline_lineage_score','optimized_fate_score','optimized_lineage_score');


%%


%plot results

%extra
[celltypecolors,~,~] = colorArrays(embeddingStruct,celltypes,lineages);

figure

subplot(1,3,1)
scatter(baseline_reduced(:,1),baseline_reduced(:,2),20,celltypecolors)
xlabel('UMAP-1')
ylabel('UMAP-2')
title('Baseline')
pbaspect([1 1 1])

subplot(1,3,2)
scatter(optimized_reduced(:,1),optimized_reduced(:,2),20,celltypecolors)
xlabel('UMAP-1')
ylabel('UMAP-2')
title('Optimized')
pbaspect([1 1 1])

baseline = [baseline_fate_score,baseline_lineage_score];
optimized = [optimized_fate_score,optimized_lineage_score];

both_scores = [baseline',optimized']';

names_chars = {'baseline','optimized'};
names_cat = categorical(names_chars);
names_cat = reordercats(names_cat,names_chars);

subplot(1,3,3)
bar(names_cat,both_scores);
title('fate score (b), lineage score (r)')

sgtitle('optimization june 29 true greedy (update parameter space) for lineage continuity')



%%
function [fate_overall, lineage_overall] = computeScore3(adjusted_percellfeature,adjusted_embeddingStruct,embryo,number_neighbors,time_window,condensed_percellfeature,condensed_embeddingStruct,starttime,sampling,endtime)

    %compute scores where lineage has a time window (on per cell per time) and fate is on the
    %condensed dataset
    
    
    %compute time window for lineage (taken from computefeatures4--see for
    %more details)
    
    %value is number of minutes forward and back to look

    %   parameter value t becomes fraction of time from 4 cell to twitch ~~ 400 mins
    %   twitch may not be measurable in dataset, so approx using 350cell tpt

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

    %approximate twitch timepoint
    %   p = approximate ratio of time from 4cell to twitch over 4cell to 350 cell = (475-75)/(345-75)  (
    %   from https://www.wormatlas.org/hermaphrodite/introduction/mainframe.html)
    p = (475-75)/(345-75);
    approx_twitch_time = p*(threefiftycelltime - fourcelltime) + fourcelltime;

    %convert minutes to corresponding numbers of timesteps in embryo
    %dataset
    time_window = ceil(time_window/400 * (approx_twitch_time-fourcelltime));
   
   
    %LINEAGE CONTINUITY
    lineage_percellscores = [];
    
    for j = 1:length(adjusted_percellfeature)


        %score is  percentage of nns that are same cell, progenitor, or
        %daughter

        %find neighbors in feature space
        [~,inds] = pdist2(adjusted_percellfeature,adjusted_percellfeature(j,:),'euclidean','Smallest',number_neighbors+1);
        inds = inds(2:end);
        
        
        %find the embryo coords of cells that are in continuous lineage
        %w/ cell j
        
        continuous_cells = []; %embryo coords, first column embryo time, second column cell index
        
        
        start_time = adjusted_embeddingStruct(j).captureFrame;
        start_ind = adjusted_embeddingStruct(j).captureIndex;

        
        %obtain past cells  
        [~,~,~,~,all_prev_cell_times,all_prev_cell_inds] = backTraverse_keepall(start_time,start_ind,'timesteps',time_window,embryo);%change to an amount of time
        continuous_cells = [continuous_cells;all_prev_cell_times,all_prev_cell_inds];
        
        
        %same process; make a new function that returns all future cell
        %times, inds, map them to coords in feature space, tally
        
        %obtain future cells
        [~,~,all_later_cell_times,all_later_cell_inds] = forwardTraverse_keepall(start_time,start_ind,'timesteps',time_window,embryo);%change to an amount of time
        continuous_cells = [continuous_cells;all_later_cell_times,all_later_cell_inds];
        
        
        %convert included cells to embeddingStruct (featurespace) indices
        
        [embdatToEmbeddingStruct,~] = mappingMaker(embryo,starttime,sampling,endtime);
        
        continuous_cells_featureCoords = zeros(size(continuous_cells,1),1);
        for l = 1:size(continuous_cells,1)
            
            continuous_cells_featureCoords(l) = embdatToEmbeddingStruct(continuous_cells(l,1)).map(continuous_cells(l,2));

        end
        
        
        %tally how many neighbors are 'continuous' (in included cells)
        num_continuous = 0;
        for k = 1:length(inds)
   
            if ismember(inds(k),continuous_cells_featureCoords)
                num_continuous = num_continuous+1;
            end

        end

        %find percent continuous, add to lineage_percellscores
        percent_continuous = num_continuous/number_neighbors;
        lineage_percellscores = [lineage_percellscores;percent_continuous];

    end
    
    
    
    
    %TERMINAL FATE SIMILARITY
    
    fate_percellscores = [];
    for j = 1:length(condensed_percellfeature)
        
        %find neighbors in feature space
        [~,condensed_inds] = pdist2(condensed_percellfeature,condensed_percellfeature(j,:),'euclidean','Smallest',number_neighbors+1);
        condensed_inds = condensed_inds(2:end);

        
        if condensed_embeddingStruct(j).celltype ~= "OTHER"
            %score is  percentage of nns that are same type as cell in
            %question

            %tally how many neighbors are of same type
            num_same_type = 0;
            for k = 1:number_neighbors

                if condensed_embeddingStruct(condensed_inds(k)).celltype == condensed_embeddingStruct(j).celltype
                    num_same_type = num_same_type+1;
                end

            end
            %percent of neighbors that are the same type
            percent_same_type = num_same_type/number_neighbors;

            %log score
            fate_percellscores = [fate_percellscores;percent_same_type];

        end
        
    end
    
    
    %average (overall) scores
    fate_overall = mean(fate_percellscores);
    lineage_overall = mean(lineage_percellscores);

end

function [condensed_feature,condensed_embeddingStruct] = condenseToOneEntryPerCell(feature,embeddingStruct,mode)
    %take in a feature space/embeddingStruct where each cell is
    %represented and multiple times and return versions of both where each
    %cell is represented only once
    
    
    condensed_feature = [];
	condensed_embeddingStruct = embeddingStruct(1);%seed with first entry to avoid nonlike struct error
    c = 1;
    names_list = ["seed"];
    for i = 1:size(feature,1)
        names_list(i) = convertCharsToStrings(embeddingStruct(i).name);
    end
    queue = [1:size(feature,1)];
    
    while ~isempty(queue)
        %get name of first cell in there
        top_name = embeddingStruct(queue(1)).name;
        
        
        %get inds of all entries of that cell
        same_cell_inds = find(names_list == top_name);
        
        %make embStruct with all entries for that particular cell
        same_cellStruct = embeddingStruct(same_cell_inds);
        
        switch mode
            case 'average_of_all_entries' %conceptually doesn't make sense lol
                
                %update condensed_feature and embeddingStruct
                feat_entry = mean(feature(same_cell_inds,:),1);
                condensed_feature(c,:) = feat_entry;
                
                condensed_embeddingStruct(c) = embeddingStruct(queue(1));
                c = c+1;
                
            case 'median_timepoint_entry'
                times = [embeddingStruct(same_cell_inds).captureFrame];%correspond to same_cell_inds
                
                %find the 'median' timept entry
                sorted_times = sort(times);
                med_tpt = sorted_times(ceil(length(sorted_times)/2));%ceil ensures that it works even w/one time
                entry_to_use = same_cell_inds(find(times==med_tpt));
                
                %update condensed_feature and embeddingStruct
                feat_entry = feature(entry_to_use,:);
                condensed_feature(c,:) = feat_entry;
                condensed_embeddingStruct(c) = embeddingStruct(entry_to_use);
                c = c+1;
        end
        %remove indices from queue
        queue = setdiff(queue,same_cell_inds);
    end
    
    %remove now meaningless fields from embeddingStruct
    
    
    
    %condensed_embeddingStruct = rmfield(condensed_embeddingStruct,{'captureFrame','captureIndex'});
    
    
    
    %note: numCellsAtCapture and realSpaceCoord arent fully meaningful or
    %accurate if 'average_of_all_entries' option is used, but may be a
    %rough reference
    

    
    %{
    queue of indices to visit = all indices
    
    while queue not empty
        see what cell name it is
        find all entries with same cell name
        organize them by time
        either take avg median entry
        append to vars condensed feature, condensed_embeddingStruct
        remove indices corresponding to that cell name from queue
        
    
        reason through whether indices (from find) and queue will work
            
    %}
  
end 

function [embdatToEmbeddingStruct,embeddingStructToEmbdat] = mappingMaker(embdat,starttime,sampling,endtime)
%Makes struct mappings between embdat format data and embeddingStruct-data 

embdatToEmbeddingStruct = struct();
embeddingStructToEmbdat = struct();


times = [starttime:sampling:endtime];

c = 1;
for time_index  = 1:length(times)
    
    t = times(time_index);
    
   	map_temp = zeros(size(embdat(t).finalpoints,1),1);
    for i = 1:size(embdat(t).finalpoints,1)
        
        map_temp(i) = c;
        
        %embeddingStruct->embdat
        embeddingStructToEmbdat(c).embryoTime = t;
        embeddingStructToEmbdat(c).embryoInd = i;
 
        %increment
        c = c+1;
    end
    %embdat->embeddingStruct
    embdatToEmbeddingStruct(t).map = map_temp;
    
    
end

%populate empty embdatToEmbeddingStruct entries with arrays so that they
%can be indexed into even if they arent present in the embeddingStruct

times_notIncluded = setdiff([1:length(embdat)],times);

for i = 1:length(times_notIncluded)
    
    t = times_notIncluded(i);
    embdatToEmbeddingStruct(t).map = zeros(size(embdat(t).finalpoints,1),1);
    
end

end

function [celltypecolors,lineagecolors,timecolors] = colorArrays(embeddingStruct,celltypes,lineages)

numPointsInEmbedding = size(embeddingStruct,2);
endtime = embeddingStruct(end).captureFrame;


%color scheme

%temporal evolution of embryo
timecolors=zeros(numPointsInEmbedding,3);
for i = 1:numPointsInEmbedding
    hsvcolor = [embeddingStruct(i).captureFrame/endtime,1,1];
    timecolors(i,:)=hsv2rgb(hsvcolor);
end

%cell type 

%%%*** remove this if statement instead make the option whether
%or not you want to have a 3d plot too

celltypecategories = [[0,1,0];[0,0,1];[1,0,0];[1,0,1];[0,1,1]];

%default color gray
celltypecolors = repmat([.5,.5,.5],numPointsInEmbedding,1);

for i = 1:numPointsInEmbedding
    for j = 1:length(celltypes)
        if embeddingStruct(i).celltype == celltypes(j)
            celltypecolors(i,:) = celltypecategories(j,:);
        end
    end

end

%lineage 
%note: must be same length as lineages above
lineagecategories = [[1,0,0];[0,1,0];[0,0,1];[1,0,1];[1,.5,0];[0,1,1];[.5,0,1]];

%default color gray
lineagecolors = repmat([.25,.25,.25],numPointsInEmbedding,1);

for i = 1:numPointsInEmbedding
    for j = 1:length(lineages)
        if embeddingStruct(i).lineage == lineages(j)
            lineagecolors(i,:) = lineagecategories(j,:);
        end
    end
end

end

