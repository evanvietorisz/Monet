%tuning hyperparameters jan 14

%for each case in cases, gives score = percent of 20 nns in feature space
%that are same type as a cell in question, averaged over all terminal cells

%requires preloading
%   embdat_stabilized
%   embeddingStruct
%   starttime
%   endtime
%   sampling


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


cases = load('optimization_cases_mar13_1.mat');
cases = cases.optimization_cases_mar13_1;

save_name = 'optimization_cases_mar13_1_WITHRESULTS.mat';

times = [starttime:sampling:endtime];
normalize_bool = true;
input_xy_res = .16;
    
number_neighbors = 20; %used to compute score (in feature space)



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
    
        
%compute scores for all cases     
for i = 1:size(cases,1)
    
    %specify params
    feats = cases{i,1};
    param_values = cases{i,2};
    
    try
        %extract featurevector
        [percellfeature,normed_percellfeature,~] = computeFeatures4(embryo,times,feats,param_values,normalize_bool,input_xy_res);
        
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
    
    
    %compute scores
    
    %terminal fate--old
    %{
    %terminal fate
    fate_percellscores = [];
    
    for j = 1:length(embeddingStruct)
        
        if embeddingStruct(j).celltype ~= "OTHER"
            %score is  percentage of nns that are same type as cell in
            %question
            
            %find neighbors
            [~,inds] = pdist2(percellfeature,percellfeature(j,:),'euclidean','Smallest',number_neighbors+1);
            inds = inds(2:end);
            
            %tally how many neighbors are of same type
            num_same_type = 0;
            for k = 1:number_neighbors
                
                if embeddingStruct(inds(k)).celltype == embeddingStruct(j).celltype
                    num_same_type = num_same_type+1;
                end
                
            end
            %percent of neighbors that are the same type
            percent_same_type = num_same_type/number_neighbors;
            
            %log score
            fate_percellscores = [fate_percellscores;percent_same_type];
            
        end
        
    end
    %}
    
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

save(save_name,'cases_output')















     