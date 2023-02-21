% Jul 2--a knn classifier to predict cell type in one embryo embedding from
% the embedding of another

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
          
          
%optimizations 7/2
jul2_20nns_fate = {"radius_boundingsphere_celldensity",15;
    "radius_boundingsphere_residual",   10;
    "radius_boundingsphere_nucsize",    25;
    "k_intercelldistance",              5;
    "degree_cousins",                   1;
    "timewindow_migrationpaths",        .15;
    "timewindow_migrationcongruency",   .0375;
    "k_migrationcongruency",            5};

jul2_20nns_lineage = {"radius_boundingsphere_celldensity",25;
    "radius_boundingsphere_residual",   10;
    "radius_boundingsphere_nucsize",    5;
    "k_intercelldistance",              5;
    "degree_cousins",                   1;
    "timewindow_migrationpaths",        .15;
    "timewindow_migrationcongruency",   .0625;
    "k_migrationcongruency",            5};


jul2_5nns_fate = {"radius_boundingsphere_celldensity",25;
    "radius_boundingsphere_residual",   10;
    "radius_boundingsphere_nucsize",    25;
    "k_intercelldistance",              5;
    "degree_cousins",                   1;
    "timewindow_migrationpaths",        .1250;
    "timewindow_migrationcongruency",   .0625;
    "k_migrationcongruency",            15};


jul2_5nns_lineage = {"radius_boundingsphere_celldensity",25;
    "radius_boundingsphere_residual",   10;
    "radius_boundingsphere_nucsize",    20;
    "k_intercelldistance",              5;
    "degree_cousins",                   1;
    "timewindow_migrationpaths",        .0625;
    "timewindow_migrationcongruency",   .0750;
    "k_migrationcongruency",            15};

%%

%knn classifier trials jul 9

%{
for 20nn fate, lineage
    full D, 15 D
        5 nns, 20 nns
            swap which is training
                
%}

%load embs
emb1 = load('WT_1.mat');
emb2 = load('WT_2.mat');

%initialize
trialnum = 1;
scores = [];
descriptions = ["seed"];


%20 nn fate
params = jul2_20nns_fate;

[~,normed_percellfeature1,~] = computeFeatures4(emb1.embdat_stabilized,[emb1.starttime:emb1.sampling:emb1.endtime],feats,params,true,emb1.input_xy_res);
[~,normed_percellfeature2,~] = computeFeatures4(emb2.embdat_stabilized,[emb2.starttime:emb2.sampling:emb2.endtime],feats,params,true,emb2.input_xy_res);

    %full D
    lower_D = size(normed_percellfeature1,2);

        %use 5 nns
        number_neighbors = 5;
        
            %emb 1 is training
            which_emb = 1;

            [tissueIDscore, cellIDscore] = knnClassification(emb1,emb2,normed_percellfeature1,normed_percellfeature2,lower_D,number_neighbors,which_emb);
            trial_description = ['trial ',num2str(trialnum),': params = jul2_20nns_fate, dimension = ',num2str(lower_D),' (full), k = ',num2str(number_neighbors),', training emb is emb',num2str(which_emb),' | tissueIDscore = ',num2str(tissueIDscore),', cellIDscore = ',num2str(cellIDscore)];
            disp(trial_description)
            
            scores(c,:) = [tissueIDscore,cellIDscore];
            descriptions(c) = convertCharsToStrings(trial_description);
            trialnum = trialnum + 1;
            
            
            %emb 2 is training
            which_emb = 2;

            [tissueIDscore, cellIDscore] = knnClassification(emb1,emb2,normed_percellfeature1,normed_percellfeature2,lower_D,number_neighbors,which_emb);
            trial_description = ['trial ',num2str(trialnum),': params = jul2_20nns_fate, dimension = ',num2str(lower_D),' (full), k = ',num2str(number_neighbors),', training emb is emb',num2str(which_emb),' | tissueIDscore = ',num2str(tissueIDscore),', cellIDscore = ',num2str(cellIDscore)];
            disp(trial_description)
            
            scores(c,:) = [tissueIDscore,cellIDscore];
            descriptions(c) = convertCharsToStrings(trial_description);
            trialnum = trialnum + 1;
            
        %use 20 nns
        number_neighbors = 20;
        
            %emb 1 is training
            which_emb = 1;

            [tissueIDscore, cellIDscore] = knnClassification(emb1,emb2,normed_percellfeature1,normed_percellfeature2,lower_D,number_neighbors,which_emb);
            trial_description = ['trial ',num2str(trialnum),': params = jul2_20nns_fate, dimension = ',num2str(lower_D),' (full), k = ',num2str(number_neighbors),', training emb is emb',num2str(which_emb),' | tissueIDscore = ',num2str(tissueIDscore),', cellIDscore = ',num2str(cellIDscore)];
            disp(trial_description)
            
            scores(c,:) = [tissueIDscore,cellIDscore];
            descriptions(c) = convertCharsToStrings(trial_description);
            trialnum = trialnum + 1;
            
            
            %emb 2 is training
            which_emb = 2;

            [tissueIDscore, cellIDscore] = knnClassification(emb1,emb2,normed_percellfeature1,normed_percellfeature2,lower_D,number_neighbors,which_emb);
            trial_description = ['trial ',num2str(trialnum),': params = jul2_20nns_fate, dimension = ',num2str(lower_D),' (full), k = ',num2str(number_neighbors),', training emb is emb',num2str(which_emb),' | tissueIDscore = ',num2str(tissueIDscore),', cellIDscore = ',num2str(cellIDscore)];
            disp(trial_description)
            
            scores(c,:) = [tissueIDscore,cellIDscore];
            descriptions(c) = convertCharsToStrings(trial_description);
            trialnum = trialnum + 1;
            
    %15 D
    lower_D = 15;
    
        %use 5 nns
        number_neighbors = 5;
        
            %emb 1 is training
            which_emb = 1;

            [tissueIDscore, cellIDscore] = knnClassification(emb1,emb2,normed_percellfeature1,normed_percellfeature2,lower_D,number_neighbors,which_emb);
            trial_description = ['trial ',num2str(trialnum),': params = jul2_20nns_fate, dimension = ',num2str(lower_D),', k = ',num2str(number_neighbors),', training emb is emb',num2str(which_emb),' | tissueIDscore = ',num2str(tissueIDscore),', cellIDscore = ',num2str(cellIDscore)];
            disp(trial_description)
            
            scores(c,:) = [tissueIDscore,cellIDscore];
            descriptions(c) = convertCharsToStrings(trial_description);
            trialnum = trialnum + 1;
            
            
            %emb 2 is training
            which_emb = 2;

            [tissueIDscore, cellIDscore] = knnClassification(emb1,emb2,normed_percellfeature1,normed_percellfeature2,lower_D,number_neighbors,which_emb);
            trial_description = ['trial ',num2str(trialnum),': params = jul2_20nns_fate, dimension = ',num2str(lower_D),', k = ',num2str(number_neighbors),', training emb is emb',num2str(which_emb),' | tissueIDscore = ',num2str(tissueIDscore),', cellIDscore = ',num2str(cellIDscore)];
            disp(trial_description)
            
            scores(c,:) = [tissueIDscore,cellIDscore];
            descriptions(c) = convertCharsToStrings(trial_description);
            trialnum = trialnum + 1;
            
        %use 20 nns
        number_neighbors = 20;
        
            %emb 1 is training
            which_emb = 1;

            [tissueIDscore, cellIDscore] = knnClassification(emb1,emb2,normed_percellfeature1,normed_percellfeature2,lower_D,number_neighbors,which_emb);
            trial_description = ['trial ',num2str(trialnum),': params = jul2_20nns_fate, dimension = ',num2str(lower_D),', k = ',num2str(number_neighbors),', training emb is emb',num2str(which_emb),' | tissueIDscore = ',num2str(tissueIDscore),', cellIDscore = ',num2str(cellIDscore)];
            disp(trial_description)
            
            scores(c,:) = [tissueIDscore,cellIDscore];
            descriptions(c) = convertCharsToStrings(trial_description);
            trialnum = trialnum + 1;
            
            
            %emb 2 is training
            which_emb = 2;

            [tissueIDscore, cellIDscore] = knnClassification(emb1,emb2,normed_percellfeature1,normed_percellfeature2,lower_D,number_neighbors,which_emb);
            trial_description = ['trial ',num2str(trialnum),': params = jul2_20nns_fate, dimension = ',num2str(lower_D),', k = ',num2str(number_neighbors),', training emb is emb',num2str(which_emb),' | tissueIDscore = ',num2str(tissueIDscore),', cellIDscore = ',num2str(cellIDscore)];
            disp(trial_description)
            
            scores(c,:) = [tissueIDscore,cellIDscore];
            descriptions(c) = convertCharsToStrings(trial_description);
            trialnum = trialnum + 1;
            
            
            
            
%20 nn fate
params = jul2_20nns_lineage;

[~,normed_percellfeature1,~] = computeFeatures4(emb1.embdat_stabilized,[emb1.starttime:emb1.sampling:emb1.endtime],feats,params,true,emb1.input_xy_res);
[~,normed_percellfeature2,~] = computeFeatures4(emb2.embdat_stabilized,[emb2.starttime:emb2.sampling:emb2.endtime],feats,params,true,emb2.input_xy_res);

    %full D
    lower_D = size(normed_percellfeature1,2);

        %use 5 nns
        number_neighbors = 5;
        
            %emb 1 is training
            which_emb = 1;

            [tissueIDscore, cellIDscore] = knnClassification(emb1,emb2,normed_percellfeature1,normed_percellfeature2,lower_D,number_neighbors,which_emb);
            trial_description = ['trial ',num2str(trialnum),': params = jul2_20nns_lineage, dimension = ',num2str(lower_D),' (full), k = ',num2str(number_neighbors),', training emb is emb',num2str(which_emb),' | tissueIDscore = ',num2str(tissueIDscore),', cellIDscore = ',num2str(cellIDscore)];
            disp(trial_description)
            
            scores(c,:) = [tissueIDscore,cellIDscore];
            descriptions(c) = convertCharsToStrings(trial_description);
            trialnum = trialnum + 1;
            
            
            %emb 2 is training
            which_emb = 2;

            [tissueIDscore, cellIDscore] = knnClassification(emb1,emb2,normed_percellfeature1,normed_percellfeature2,lower_D,number_neighbors,which_emb);
            trial_description = ['trial ',num2str(trialnum),': params = jul2_20nns_lineage, dimension = ',num2str(lower_D),' (full), k = ',num2str(number_neighbors),', training emb is emb',num2str(which_emb),' | tissueIDscore = ',num2str(tissueIDscore),', cellIDscore = ',num2str(cellIDscore)];
            disp(trial_description)
            
            scores(c,:) = [tissueIDscore,cellIDscore];
            descriptions(c) = convertCharsToStrings(trial_description);
            trialnum = trialnum + 1;
            
        %use 20 nns
        number_neighbors = 20;
        
            %emb 1 is training
            which_emb = 1;

            [tissueIDscore, cellIDscore] = knnClassification(emb1,emb2,normed_percellfeature1,normed_percellfeature2,lower_D,number_neighbors,which_emb);
            trial_description = ['trial ',num2str(trialnum),': params = jul2_20nns_lineage, dimension = ',num2str(lower_D),' (full), k = ',num2str(number_neighbors),', training emb is emb',num2str(which_emb),' | tissueIDscore = ',num2str(tissueIDscore),', cellIDscore = ',num2str(cellIDscore)];
            disp(trial_description)
            
            scores(c,:) = [tissueIDscore,cellIDscore];
            descriptions(c) = convertCharsToStrings(trial_description);
            trialnum = trialnum + 1;
            
            
            %emb 2 is training
            which_emb = 2;

            [tissueIDscore, cellIDscore] = knnClassification(emb1,emb2,normed_percellfeature1,normed_percellfeature2,lower_D,number_neighbors,which_emb);
            trial_description = ['trial ',num2str(trialnum),': params = jul2_20nns_lineage, dimension = ',num2str(lower_D),' (full), k = ',num2str(number_neighbors),', training emb is emb',num2str(which_emb),' | tissueIDscore = ',num2str(tissueIDscore),', cellIDscore = ',num2str(cellIDscore)];
            disp(trial_description)
            
            scores(c,:) = [tissueIDscore,cellIDscore];
            descriptions(c) = convertCharsToStrings(trial_description);
            trialnum = trialnum + 1;
            
    %15 D
    lower_D = 15;
    
        %use 5 nns
        number_neighbors = 5;
        
            %emb 1 is training
            which_emb = 1;

            [tissueIDscore, cellIDscore] = knnClassification(emb1,emb2,normed_percellfeature1,normed_percellfeature2,lower_D,number_neighbors,which_emb);
            trial_description = ['trial ',num2str(trialnum),': params = jul2_20nns_lineage, dimension = ',num2str(lower_D),', k = ',num2str(number_neighbors),', training emb is emb',num2str(which_emb),' | tissueIDscore = ',num2str(tissueIDscore),', cellIDscore = ',num2str(cellIDscore)];
            disp(trial_description)
            
            scores(c,:) = [tissueIDscore,cellIDscore];
            descriptions(c) = convertCharsToStrings(trial_description);
            trialnum = trialnum + 1;
            
            
            %emb 2 is training
            which_emb = 2;

            [tissueIDscore, cellIDscore] = knnClassification(emb1,emb2,normed_percellfeature1,normed_percellfeature2,lower_D,number_neighbors,which_emb);
            trial_description = ['trial ',num2str(trialnum),': params = jul2_20nns_lineage, dimension = ',num2str(lower_D),', k = ',num2str(number_neighbors),', training emb is emb',num2str(which_emb),' | tissueIDscore = ',num2str(tissueIDscore),', cellIDscore = ',num2str(cellIDscore)];
            disp(trial_description)
            
            scores(c,:) = [tissueIDscore,cellIDscore];
            descriptions(c) = convertCharsToStrings(trial_description);
            trialnum = trialnum + 1;
            
        %use 20 nns
        number_neighbors = 20;
        
            %emb 1 is training
            which_emb = 1;

            [tissueIDscore, cellIDscore] = knnClassification(emb1,emb2,normed_percellfeature1,normed_percellfeature2,lower_D,number_neighbors,which_emb);
            trial_description = ['trial ',num2str(trialnum),': params = jul2_20nns_lineage, dimension = ',num2str(lower_D),', k = ',num2str(number_neighbors),', training emb is emb',num2str(which_emb),' | tissueIDscore = ',num2str(tissueIDscore),', cellIDscore = ',num2str(cellIDscore)];
            disp(trial_description)
            
            scores(c,:) = [tissueIDscore,cellIDscore];
            descriptions(c) = convertCharsToStrings(trial_description);
            trialnum = trialnum + 1;
            
            
            %emb 2 is training
            which_emb = 2;

            [tissueIDscore, cellIDscore] = knnClassification(emb1,emb2,normed_percellfeature1,normed_percellfeature2,lower_D,number_neighbors,which_emb);
            trial_description = ['trial ',num2str(trialnum),': params = jul2_20nns_lineage, dimension = ',num2str(lower_D),', k = ',num2str(number_neighbors),', training emb is emb',num2str(which_emb),' | tissueIDscore = ',num2str(tissueIDscore),', cellIDscore = ',num2str(cellIDscore)];
            disp(trial_description)
            
            scores(c,:) = [tissueIDscore,cellIDscore];
            descriptions(c) = convertCharsToStrings(trial_description);
            trialnum = trialnum + 1;


%save

save('knn_classifier_jul9.mat','scores','descriptions')

            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return
%%

%create the feature representations you want to look at

params = jul2_5nns_lineage;

lower_D = 15; %dimension of space to find knns
neighborparam = 30; %umap parameter

%both loaded from 1 (min) cells to 571,573 cells, respectively
%   the fact that neither dataset is truncated means the two separately
%   normalized datasets should still line up
emb1 = load('WT_1.mat');
emb2 = load('WT_2.mat');

%compute features
[~,normed_percellfeature1,~] = computeFeatures4(emb1.embdat_stabilized,[emb1.starttime:emb1.sampling:emb1.endtime],feats,params,true,emb1.input_xy_res);
[~,normed_percellfeature2,~] = computeFeatures4(emb2.embdat_stabilized,[emb2.starttime:emb2.sampling:emb2.endtime],feats,params,true,emb2.input_xy_res);




%reduce dimension of feature space
if lower_D == size(normed_percellfeature1,2)
    both_reduced = [normed_percellfeature1;normed_percellfeature2];
else
    both_reduced = run_umap([normed_percellfeature1;normed_percellfeature2],'n_components',lower_D,'n_neighbors',neighborparam);
end

reduced1 = both_reduced(1:length(normed_percellfeature1),:);
reduced2 = both_reduced(length(normed_percellfeature1)+1:end,:);


%%

%%% parameters %%%

which_emb = 2;%which emb is the training? 1 or 2

number_neighbors = 5; %number of neighbors to actually look at

%%% backend %%%

%assign training,validation
switch which_emb
    case 1
        training = reduced1;
        validation = reduced2;
        
        training_emb = emb1;
        validation_emb = emb2;
        
    case 2
        training = reduced2;
        validation = reduced1;
        
        training_emb = emb2;
        validation_emb = emb1;  
end

%initialize
overall_name_score = 0;
term_fate_score = 0;
num_term_cells = 0;
%%
%make lists of types,names
training_celltypes = repmat([" "],size(training,1),1);
training_cellnames = repmat([" "],size(training,1),1);

for i = 1:size(training,1)
    training_cellnames(i) = convertCharsToStrings(training_emb.embeddingStruct(i).name);
    training_celltypes(i) = convertCharsToStrings(training_emb.embeddingStruct(i).celltype);
end


%%% score %%%

for i = 1:size(validation,1)
    
    %compute neighbors
    [~,inds] = pdist2(training,validation(i,:),'euclidean','Smallest',number_neighbors+1);
    inds = inds(2:end);
    
    %obtain most frequent cell type
    most_frequent_celltype = nn_mode(training_celltypes,inds);
    
    %obtain most frequent cell name
    most_frequent_cellname = nn_mode(training_cellnames,inds);
    
    
    %get cell i's info
    own_name = validation_emb.embeddingStruct(i).name;
    own_type = validation_emb.embeddingStruct(i).celltype;
    
    %log scores
    if own_name == most_frequent_cellname
        overall_name_score = overall_name_score+1;
    end
    
    if own_type ~= 'OTHER'
        num_term_cells = num_term_cells+1;
        if own_type == most_frequent_celltype
            term_fate_score = term_fate_score+1;
        end
    end

end

%normalize to percent
overall_name_score = overall_name_score/size(validation,1); %percent of all cells for which mode of nn cellnames predicts own name
term_fate_score = term_fate_score/num_term_cells; %percent of all terminal cells "" "" celltypes "" "" own type

disp(['terminal fate score: ',num2str(term_fate_score)])
disp(['overall name score: ',num2str(overall_name_score)])
%%

function [mode_entry] = nn_mode(all_entries,inds)

    %gets mode entry of a list
    %assume inds are nearest neighbor inds in increasing order of distance
    %obtains single mode by discarding farthest neighbor until there is
    %only one mode

    %obtain most frequent cell name
    while true
        entries_to_look_at = all_entries(inds);
        [~,~,number_rep] = unique(entries_to_look_at);
        unique_entries = unique(number_rep);
        count = histc(number_rep,unique_entries);
        modes = unique_entries(count==max(count));
        
        if length(modes) == 1
            break
        end
        
        inds = inds(1:end-1);%look at fewer and fewer neighbors until there is one mode

    end
    mode_entry = entries_to_look_at(number_rep == modes);
    mode_entry = mode_entry(1);%they're all the same, so just take first

end


function [tissueIDscore, cellIDscore] = knnClassification(emb1,emb2,normed_percellfeature1,normed_percellfeature2,lower_D,number_neighbors,which_emb)

%reduce dimension of feature space
if lower_D == size(normed_percellfeature1,2)
    both_reduced = [normed_percellfeature1;normed_percellfeature2];
else
    both_reduced = run_umap([normed_percellfeature1;normed_percellfeature2],'n_components',lower_D,'n_neighbors',neighborparam);
end

reduced1 = both_reduced(1:length(normed_percellfeature1),:);
reduced2 = both_reduced(length(normed_percellfeature1)+1:end,:);


%assign training,validation
switch which_emb
    case 1
        training = reduced1;
        validation = reduced2;
        
        training_emb = emb1;
        validation_emb = emb2;
        
    case 2
        training = reduced2;
        validation = reduced1;
        
        training_emb = emb2;
        validation_emb = emb1;  
end

%make lists of types,names
training_celltypes = repmat([" "],size(training,1),1);
training_cellnames = repmat([" "],size(training,1),1);

%initialize
overall_name_score = 0;
term_fate_score = 0;
num_term_cells = 0;


%score
for i = 1:size(validation,1)
    
    %compute neighbors
    [~,inds] = pdist2(training,validation(i,:),'euclidean','Smallest',number_neighbors+1);
    inds = inds(2:end);
    
    %obtain most frequent cell type
    most_frequent_celltype = nn_mode(training_celltypes,inds);
    
    %obtain most frequent cell name
    most_frequent_cellname = nn_mode(training_cellnames,inds);
    
    
    %get cell i's info
    own_name = validation_emb.embeddingStruct(i).name;
    own_type = validation_emb.embeddingStruct(i).celltype;
    
    %log scores
    if own_name == most_frequent_cellname
        overall_name_score = overall_name_score+1;
    end
    
    if own_type ~= 'OTHER'
        num_term_cells = num_term_cells+1;
        if own_type == most_frequent_celltype
            term_fate_score = term_fate_score+1;
        end
    end

end

%normalize to percents
cellIDscore = overall_name_score/size(validation,1); %percent of all cells for which mode of nn cellnames predicts own name
tissueIDscore = term_fate_score/num_term_cells; %percent of all terminal cells "" "" celltypes "" "" own type

end

