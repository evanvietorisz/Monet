% a knn classifier to predict cell type in one embryo embedding from
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

%sequential trials using three embryos

%parameters

%dimensions into which space is reduced before classification
num_dimensions = [55,25,15,7,4,3,2];%largest must be <= horizontal dimension of smallest space (55, usually)

%number of neighbors for classification
nns_for_classification = [40,20,5,3,1];

%which embryos will be used as training *within a given pair), i.e. if you
%want only the first emb (emb1 below) in a given pair to be used as
%training, make this 1
training_emb_options = [1,2];

params = jul2_20nns_lineage;
params_name = 'jul2_20nns_lineage';

%names to save results to
trial_names = ["knn_classifier_dispim_embs_1_and_2_aug23.mat";
               "knn_classifier_dispim_embs_1_and_3_aug23.mat";
               "knn_classifier_dispim_embs_2_and_3_aug23.mat"];

%load the three embryos
embyro1 = load('dispim_emb1');
embryo2 = load('dispim_emb2');
embryo3 = load('dispim_emb3');

%perform classification
for experiment_num = 1:3
    
    %assign the right embs
    switch experiment_num
        case 1
            emb1 = embyro1; 
            emb2 = embryo2;
        case 2
            emb1 = embyro1;
            emb2 = embryo3;
        case 3
            emb1 = embryo2;
            emb2 = embryo3;
    end
    

    percellfeature1 = emb1.percellfeature;
    normed_percellfeature1 = emb1.normed_percellfeature;

    percellfeature2 = emb2.percellfeature;
    normed_percellfeature2 = emb2.normed_percellfeature;


    trial_name = trial_names(experiment_num);
    



    %initialize
    trialnum = 1;
    scores = [];
    descriptions = ["seed"];
    trials_cell = cell(1);%columns: params_name,feature space dimension, k for classification,training emb number,tissueIDscore,cellIDscore

    %do trial
    for i = 1:length(num_dimensions)
        for j = 1:length(nns_for_classification)
            for k = 1:length(training_emb_options)

                lower_D = num_dimensions(i);
                number_neighbors = nns_for_classification(j);
                which_emb = training_emb_options(k);

                [tissueIDscore, cellIDscore] = knnClassification(emb1,emb2,normed_percellfeature1,normed_percellfeature2,lower_D,number_neighbors,which_emb);

                trial_description = ['trial ',num2str(trialnum),': params = ',params_name,', dimension = ',num2str(lower_D),', k = ',num2str(number_neighbors),', training emb is emb',num2str(which_emb),' | tissueIDscore = ',num2str(tissueIDscore),', cellIDscore = ',num2str(cellIDscore)];
                disp(trial_description)

                trials_cell{trialnum,1} = params_name;
                trials_cell{trialnum,2} = lower_D;
                trials_cell{trialnum,3} = number_neighbors;
                trials_cell{trialnum,4} = which_emb;
                trials_cell{trialnum,5} = tissueIDscore;
                trials_cell{trialnum,6} = cellIDscore;

                scores(trialnum,:) = [tissueIDscore,cellIDscore];
                descriptions(trialnum) = convertCharsToStrings(trial_description);
                trialnum = trialnum + 1;

            end 
        end
    end

    save(trial_name,'scores','descriptions','trials_cell')


end







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
    both_reduced = run_umap([normed_percellfeature1;normed_percellfeature2],'n_components',lower_D,'n_neighbors',30);
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

%make lists of types,names to index into w/ inds
training_celltypes = repmat([" "],size(training,1),1);
training_cellnames = repmat([" "],size(training,1),1);

for i = 1:size(training,1)
    training_celltypes(i) = convertStringsToChars(training_emb.embeddingStruct(i).celltype);
    training_cellnames(i) = convertStringsToChars(training_emb.embeddingStruct(i).name);
end




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

function [normed_featurevector] = quick_normalize(featurevector,indices_list)
    %* scales time features from min to max within the dataset, not using
    %some external landmark like 4 cells stage--this means any two feature
    %vectors being compared using this function need to span the same
    %developmental window (in terms of num cells in embryo)


    %*assumes that all feats are present and ordered as in:
    %feature list
    %{

    1     "beginning_of_current_cell_cycle"
    2     "end_of_current_cell_cycle"
    3     "duration_of_current_cell_cycle"
    4     "time_until_next_division"
    5     "time_since_last_division"
    6     "fraction_of_current_cycle_elapsed"
    7     "similarity_to_ancestors_cell_cycle"
    8     "similarity_to_cousins_cell_cycle"
    9     "first_deriv_x"
    10    "first_deriv_y"
    11    "first_deriv_z"
    12    "laplacian"
    13    "local_normal_vector"
    14    "residual"
    15    "distances_to_k_nns"
    16    "signed_distances"
    17    "vec_to_position_x_mins_ago"
    18    "vec_to_position_x_mins_in_future"
    19    "vec_to_site_of_last_division"
    20    "vec_to_site_of_next_division"
    21    "displacement_vec_over_cell_cycle"
    22    "proportion_of_present_nns_that_are_descendants_of_past_nns"
    23    "filtered_nuc_size"
    24    "number_of_cell_in_embryo"
    25    "developmental_time"

    %}


    %normalizations corresponding to feats
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

    %make copy
    normed_featurevector = featurevector;

    for i = 1:size(indices_list,1)

        %get piece of normed_featurevector
        entries = normed_featurevector(:,indices_list(i,1):indices_list(i,2));

        if ordered_norms(i) == "scale_4cell_to_end"
            %linearly scale all entries on [0,1] (using range of values in dataset, not external 4 cell stage)
            entries = rescale(entries);

        elseif ordered_norms(i) == "scale_0_to_max_num_cells"
            %linearly scale all entries on [0,1] (using range of values in dataset, not external 4 cell stage)
            entries = rescale(entries);

        elseif ordered_norms(i) == "column_zscore"
            %take zscore of column ( indices_list(i,1)==indices_list(i,2) )
            entries = zscore(entries);

        elseif ordered_norms(i) == "all_zscore"
            %take zscore of all elements
            entries = zscore(entries,0,'all');

        elseif ordered_norms(i) == "divide_avg_mag"
            %divide array elements by avg magitude of a vector in the array
            entries = avgMagnitude(entries);
        end

        %put newly normalized entries back in normed_featurevector
        normed_featurevector(:,indices_list(i,1):indices_list(i,2)) = entries;

    end

end
