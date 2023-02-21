function [colors,seedCoordsArray,seedConnectorsArray] = backpropagateLineages(seedcell,connectorsBool,embeddingStruct)
%color lineages backpropogated from user-specified terminal cells

%seedcell is cell list of cell names (char arrays) to backpropogate as distinct lineages
% ex. {'Caaaaa';'Caaaap';'Caaapa';'Caaapp'}...
%taken directly from .mat file

%connectorsBool is a bool deciding whether to have arrows connecting the
%points in the individual lineages

%embeddingStruct is embeddingStruct

%colors is color array for scatter plot
%seedCoordsArray is the positions of the all the points in the specified
%lineages

%note: will have an error assigning colors for more than 20 distinct lineages


%aux
seeds = convertCharsToStrings(seedcell);
numPointsInEmbedding = size(embeddingStruct,2);

%load list of 20 distinguishable colors from
%https://sashamaps.net/docs/resources/20-colors/
load('twentycolors.mat','twentycolors');



%create color array to be used for scatter plot
colors = repmat([.5,.5,.5],numPointsInEmbedding,1);

%collect info for connection arrows
%arrows show discrete
seedCoords = cell(length(seeds),1);

%go backwards from end
for i = numPointsInEmbedding:-1:1
    %look at each seed
    for j = 1:length(seeds)
        %if cell is seed or ancestor of seed
        if startsWith(seeds(j),string(embeddingStruct(i).name))
            %mark,redefine seed
            colors(i,:) = twentycolors(j,:);
            seedCoords{j} = [seedCoords{j};embeddingStruct(i).twoDpoint];
            seeds(j) = string(embeddingStruct(i).name);
        end
    end
end

%positions of seeds in each lineage and associated connectors
seedCoordsArray = [];
seedConnectorsArray = [];

seedConnectors = cell(size(seedCoords));
%get quiver coords to show specific lineages
for i = 1:length(seedCoords)
    %seedCoords starts terminal and goes backwards; first (terminal entry has no arrow)
    seedConnectors{i} = [0,0];
    for j = 2:size(seedCoords{i},1)
        seedConnectors{i} = [seedConnectors{i};seedCoords{i}(j-1,:)-seedCoords{i}(j,:)];
    end
end

%concatenate seed coordinates and corresponding arrows into long matrices
for i = 1:length(seedCoords)
    seedCoordsArray = [seedCoordsArray;seedCoords{i}];
    seedConnectorsArray = [seedConnectorsArray;seedConnectors{i}];
end

%if you don't want connecting
if connectorsBool == false
    seedConnectorsArray = zeros(size(seedConnectorsArray));
end


end

