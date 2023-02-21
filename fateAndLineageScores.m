function [fate_overall,lineage_overall,fate_percellscores,lineage_percellscores] = fateAndLineageScores(embryo,percellfeature,embeddingStruct,number_neighbors)
%calculates the fate clustering and lineage continuity scores--see notes
%for what each line computes
% embryo is the embdat dataset
% percellfeature is the featurevector
% embeddingStruct is the embeddingStruct (must have indices that correspond
% to percellfeature
% embeddingStruct is embeddingStruct

%number_neighbors is the number of neighbors 



%compute scores

    lineage_percellscores = [];
    fate_percellscores = [];
    
    for j = 1:length(percellfeature)
        
        %LINEAGE CONTINUITY

        %score is  percentage of nns that are same cell, progenitor, or
        %daughter
        
        %find neighbors
        [~,inds] = pdist2(percellfeature,percellfeature(j,:),'euclidean','Smallest',number_neighbors+1);
        inds = inds(2:end);

        %create list of names that qualify as same cell, progenitor, or
        %daughter
        
        %add same cell
        correct_names = [convertCharsToStrings(embeddingStruct(j).name)];
        
        start_time = embeddingStruct(j).captureFrame;
        start_ind = embeddingStruct(j).captureIndex;
        
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

            if ~isempty(find(correct_names == embeddingStruct(inds(k)).name))
                num_continuous = num_continuous+1;
            end

        end
        
        %find percent continuous, add to lineage_percellscores
        percent_continuous = num_continuous/number_neighbors;
        lineage_percellscores = [lineage_percellscores;percent_continuous];
        
        
        %TERMINAL FATE SIMILARITY
        if embeddingStruct(j).celltype ~= "OTHER"
            %score is  percentage of nns that are same type as cell in
            %question
            
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
    
    
    %output
    
    fate_overall = mean(fate_percellscores);
    lineage_overall = mean(lineage_percellscores);

end

