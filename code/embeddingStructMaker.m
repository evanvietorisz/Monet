function [embeddingStruct] = embeddingStructMaker(embdat_stabilized,starttime,sampling,endtime,lineages)
%Makes a struct of the embeddings for easy lookup

%reducedfeature is 2d reduced feature
%reducedfeature3 is 3d reduced feature

%embdat_stabilized is aligned embryo struct

%term_fates_to_keep is cell array of char vectors that are the terminal
%fates of all desired cells
%NOTE: if term_fates_to_keep is empty, then ALL cells will be kept

%starttime, sampling, endtime are start, increment, and end of feature
%extraction window
%lineages is array of strings which beginning characters of desired
%lineages
%ex. lineages = ["ABa","M","E","P"]

%load cells
load('allpharynx.mat','allpharynx');
load('allneuron.mat','allneuron');
load('allbodywallmuscle.mat','allbodywallmuscle');
load('allhypoderm.mat','allhypoderm');
load('allgut.mat','allgut');

embeddingStruct = struct([]);


times = [starttime:sampling:endtime];

c = 1;
for time_index  = 1:length(times)
    
    t = times(time_index);
    
    for i = 1:size(embdat_stabilized(t).finalpoints,1) %length of reducedfeature

        
        
        %reduced coordinates
        %embeddingStruct(c).twoDpoint = reducedfeature(c,:);
        %embeddingStruct(c).threeDpoint = reducedfeature3(c,:);
        
        %cell name (char array)
        name = embdat_stabilized(t).names{i,1};
        embeddingStruct(c).name = name;
        
        %cell type (string)
        embeddingStruct(c).celltype = "OTHER";
        if ~isempty(find(strcmp(name,allpharynx), 1))
            embeddingStruct(c).celltype = "pharynx";
        elseif ~isempty(find(strcmp(name,allneuron), 1))
            embeddingStruct(c).celltype = "neuron";
        elseif ~isempty(find(strcmp(name,allbodywallmuscle), 1))
            embeddingStruct(c).celltype = "bodywallmuscle";
        elseif ~isempty(find(strcmp(name,allhypoderm), 1))
            embeddingStruct(c).celltype = "hypoderm";
        elseif ~isempty(find(strcmp(name,allgut), 1))
            embeddingStruct(c).celltype = "gut";
        end
        
        %lineage
        embeddingStruct(c).lineage = "OTHER";
        for j = 1:length(lineages)
            if startsWith(name,lineages(j))
                embeddingStruct(c).lineage = lineages(j);
            end
        end
  
        %number of cells in emb at time of capture (developmental time)
        embeddingStruct(c).numCellsAtCapture = size(embdat_stabilized(t).finalpoints,1);
        
        %real space point
        embeddingStruct(c).realSpaceCoord = embdat_stabilized(t).finalpoints(i,:);
        
        %frame of capture (for assessing continuity in objective function)
        embeddingStruct(c).captureFrame = t;
        embeddingStruct(c).timeIndex = time_index;
        embeddingStruct(c).captureIndex = i;
        
        
        %increment
        c = c+1;
    end
end

end

