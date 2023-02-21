function [scoreStruct2d,scoreStruct3d] = objectiveFunction2(embeddingStruct,celltypes,lineages)
%objective function to score plots
%embeddingStruct is embeddingStruct
%celltypes is array of strings which are the cell types
%lineages is array of strings which are the lineages

%{
Outputs:

celltype-,lineage-,frameNumSameInNeighborhood:
    For all non 'other' cells, the average percentage of n nearest
    neighbors (in embedding space) which are the same
    celltype,lineage,captureFrame* as it, respectively

    *in order to capture continuity, this includes cells with the same
    capture frame or the frame immediately before or after the cell in
    question

celltype-,lineage-,frameSmallestDistToOther
    For all on 'other' cells, the average distance to the nearest cell of a
    different type, normalized by the average distance to the nearest
    neighbor cell for all cells in the embedding
    
celltype-,lineage-,frameCompactness
    For each cell type, lineage, or captureFrame in the dataset, the
    standard deviation of distance from all cells of a given type to their
    centroid, normalized by the standard deviation of the distance from all
    cells in the embedding to their centroid 
%}


%indices correspond to indices in embeddingStruct
scoreStruct2d = struct;
scoreStruct3d = struct;

%%% aux quantities %%%
numPointsInEmbedding = size(embeddingStruct,2);

%number of neighbors to look at
num_neighbors = round(numPointsInEmbedding/1000);

%grab sampling
sampling = embeddingStruct(2).captureFrame - embeddingStruct(1).captureFrame;

%array of distnct frames
distinct_frames = [embeddingStruct(1).captureFrame];
for i = 2:numPointsInEmbedding
    if embeddingStruct(i).captureFrame ~= distinct_frames(end)
        distinct_frames = [distinct_frames;embeddingStruct(i).captureFrame];
    end
end

%get locations of all embeddings in array format
points2d = [];
for i = 1:numPointsInEmbedding
    points2d = [points2d;embeddingStruct(i).twoDpoint(1,:)];
end

points3d = [];
for i = 1:numPointsInEmbedding
    points3d = [points3d;embeddingStruct(i).threeDpoint(1,:)];
end

%average distance to nearest neighbor in whole embedding
dist_to_neighbor = 0;
for i = 1:numPointsInEmbedding
    pt = embeddingStruct(i).twoDpoint;
    [dist,ind] = pdist2(points2d,pt,'euclidean','Smallest',2);
    dist_to_neighbor = dist_to_neighbor + dist(2);
end
dist_to_neighbor = dist_to_neighbor/numPointsInEmbedding;





%%% compute scores %%% 

%proximity of like points: celltype, lineage, capture frame
celltypeNumSameInNeighborhood = 0;
lineageNumSameInNeighborhood = 0;
frameNumSameInNeighborhood = 0;

smallest_dist_to_other_celltype = 0;
smallest_dist_to_other_lineage = 0;
smallest_dist_to_other_frame = 0;

num_tallied = 0;
for i = 1:numPointsInEmbedding
    if embeddingStruct(i).celltype ~= "OTHER"
        num_tallied = num_tallied + 1;
        
        %2d
        %nearest n points in embedding
        [dists,inds] = pdist2(points2d,embeddingStruct(i).twoDpoint(1,:),'euclidean','Smallest',num_neighbors+1);
        inds = inds(2:end,:);
        
        celltype_num_same = 0;
        lineage_num_same = 0;
        frame_num_same = 0;
        
        %if same, tally
        for j = 1:length(inds)
            %cell type
            if embeddingStruct(inds(j)).celltype == embeddingStruct(i).celltype
                celltype_num_same = celltype_num_same + 1;
            end
            
            %lineage
            if embeddingStruct(inds(j)).lineage == embeddingStruct(i).lineage
                lineage_num_same = lineage_num_same + 1;
            end
            
            %developmental time (by frame)
            frameDifference =  embeddingStruct(inds(j)).captureFrame - embeddingStruct(i).captureFrame;
            if abs(frameDifference) <= sampling
                frame_num_same = frame_num_same + 1;
            end
        end
        
        %normalize to number of neigbors
        celltype_num_same = celltype_num_same/num_neighbors;
        lineage_num_same = lineage_num_same/num_neighbors;
        frame_num_same = frame_num_same/num_neighbors;
        
        %tally
        celltypeNumSameInNeighborhood = celltypeNumSameInNeighborhood + celltype_num_same;
        lineageNumSameInNeighborhood = lineageNumSameInNeighborhood + lineage_num_same;
        frameNumSameInNeighborhood = frameNumSameInNeighborhood + frame_num_same;
        
        %distance to nearest cell that is different from current (assuming it's in the neighborhood)
        %celltype
        for j = 1:length(inds)
            if embeddingStruct(inds(j)).celltype ~= embeddingStruct(i).celltype
                smallest_dist_to_other_celltype = smallest_dist_to_other_celltype + dists(j);
                break
            end
        end
        %lineage
        for j = 1:length(inds)
            if embeddingStruct(inds(j)).lineage ~= embeddingStruct(i).lineage
                smallest_dist_to_other_lineage = smallest_dist_to_other_lineage + dists(j);
                break
            end
        end
        %frame
        for j = 1:length(inds)
            if embeddingStruct(inds(j)).captureFrame ~= embeddingStruct(i).captureFrame
                smallest_dist_to_other_frame = smallest_dist_to_other_frame + dists(j);
                break
            end
        end
    end
end

%celltype 
%normalize to points tallied
celltypeNumSameInNeighborhood = celltypeNumSameInNeighborhood/num_tallied;
lineageNumSameInNeighborhood = lineageNumSameInNeighborhood/num_tallied;
frameNumSameInNeighborhood = frameNumSameInNeighborhood/num_tallied;

%add to output struct
scoreStruct2d.celltypeNumSameInNeighborhood = celltypeNumSameInNeighborhood;
scoreStruct2d.lineageNumSameInNeighborhood = lineageNumSameInNeighborhood;
scoreStruct2d.frameNumSameInNeighborhood = frameNumSameInNeighborhood;

%distance to nearest 'other'
%normalize by num tallied
smallest_dist_to_other_celltype = smallest_dist_to_other_celltype/num_tallied;
smallest_dist_to_other_lineage = smallest_dist_to_other_lineage/num_tallied;
smallest_dist_to_other_frame = smallest_dist_to_other_frame/num_tallied;

%normalize by average distance to any neighbor
smallest_dist_to_other_celltype = smallest_dist_to_other_celltype/dist_to_neighbor;
smallest_dist_to_other_lineage = smallest_dist_to_other_lineage/dist_to_neighbor;
smallest_dist_to_other_frame = smallest_dist_to_other_frame/dist_to_neighbor;

scoreStruct2d.celltypeSmallestDistToOther = smallest_dist_to_other_celltype;
scoreStruct2d.lineageSmallestDistToOther = smallest_dist_to_other_lineage;
scoreStruct2d.frameSmallestDistToOther = smallest_dist_to_other_frame;


%cluster compactness
%celltype
celltype_deviation = 0;
types_not_present = 1;
for j = 1:length(celltypes)
    celltype_current = celltypes(j);
    
    %grab indices,positions of cells of a given type
    celltype_indices = [];
    for i = 1:numPointsInEmbedding
        if embeddingStruct(i).celltype == celltype_current
            celltype_indices = [celltype_indices;i];
        end
    end
    celltype_points = points2d(celltype_indices,:);
    
    %prevent NaN values
    if isempty(celltype_points)
        types_not_present = types_not_present + 1;
        continue
    end
    
    %distances to centroid of celltype cluster
    centroid = [mean(celltype_points(:,1)),mean(celltype_points(:,2))];
    dists = pdist2(centroid,celltype_points);
    
    %log
    celltype_deviation = celltype_deviation + std(dists);
end
%normalize by number of organs
celltype_deviation = celltype_deviation/(length(celltypes)-types_not_present);


%normalize by same treatment on whole embedding
overall_centroid = [mean(points2d(:,1)),mean(points2d(:,2))];
overall_dists = pdist2(overall_centroid,points2d);
overall_std = std(overall_dists);
celltype_deviation = celltype_deviation/overall_std;

scoreStruct2d.celltypeCompactness = celltype_deviation;

%lineage
lineage_deviation = 0;
for j = 1:length(lineages)
    lineage_current = lineages(j);
    
    %grab indices,positions of cells of a given type
    lineage_indices = [];
    for i = 1:numPointsInEmbedding
        if embeddingStruct(i).lineage == lineage_current
            lineage_indices = [lineage_indices;i];
        end
    end
    lineage_points = points2d(lineage_indices,:);
    
    %distances to centroid of lineage cluster
    centroid = [mean(lineage_points(:,1)),mean(lineage_points(:,2))];
    dists = pdist2(centroid,lineage_points);
    
    %log
    lineage_deviation = lineage_deviation + std(dists);
end

%normalize by number of organs
lineage_deviation = lineage_deviation/length(lineages);

%normalize by same treatment on whole embedding
lineage_deviation = lineage_deviation/overall_std;

scoreStruct2d.lineageCompactness = lineage_deviation;


%frame
frame_deviation = 0;
for j = 1:length(distinct_frames)
    frame_current = distinct_frames(j);
    
    %grab indices,positions of cells of a given type
    frame_indices = [];
    for i = 1:numPointsInEmbedding
        if embeddingStruct(i).captureFrame == frame_current
            frame_indices = [frame_indices;i];
        end
    end
    frame_points = points2d(frame_indices,:);
    
    %distances to centroid of lineage cluster
    centroid = [mean(frame_points(:,1)),mean(frame_points(:,2))];
    dists = pdist2(centroid,frame_points);
    
    %log
    frame_deviation = frame_deviation + std(dists);
end

%normalize by number of organs
frame_deviation = frame_deviation/length(distinct_frames);

%normalize by same treatment on whole embedding
frame_deviation = frame_deviation/overall_std;

scoreStruct2d.frameCompactness = frame_deviation;

%%% scores 3d %%%-- terrible to repeat everything will clean up later

%proximity of like points: celltype, lineage, capture frame

%reinitialize
celltypeNumSameInNeighborhood = 0;
lineageNumSameInNeighborhood = 0;
frameNumSameInNeighborhood = 0;

smallest_dist_to_other_celltype = 0;
smallest_dist_to_other_lineage = 0;
smallest_dist_to_other_frame = 0;

num_tallied = 0;
for i = 1:numPointsInEmbedding
    if embeddingStruct(i).celltype ~= "OTHER"
        num_tallied = num_tallied + 1;
        
        %2d
        %nearest n points in embedding
        [dists,inds] = pdist2(points3d,embeddingStruct(i).threeDpoint(1,:),'euclidean','Smallest',num_neighbors+1);
        inds = inds(2:end,:);
        
        celltype_num_same = 0;
        lineage_num_same = 0;
        frame_num_same = 0;
        
        %if same, tally
        for j = 1:length(inds)
            %cell type
            if embeddingStruct(inds(j)).celltype == embeddingStruct(i).celltype
                celltype_num_same = celltype_num_same + 1;
            end
            
            %lineage
            if embeddingStruct(inds(j)).lineage == embeddingStruct(i).lineage
                lineage_num_same = lineage_num_same + 1;
            end
            
            %developmental time (by frame)
            frameDifference =  embeddingStruct(inds(j)).captureFrame - embeddingStruct(i).captureFrame;
            if abs(frameDifference) <= sampling
                frame_num_same = frame_num_same + 1;
            end
        end
        
        %normalize to number of neigbors
        celltype_num_same = celltype_num_same/num_neighbors;
        lineage_num_same = lineage_num_same/num_neighbors;
        frame_num_same = frame_num_same/num_neighbors;
        
        %tally
        celltypeNumSameInNeighborhood = celltypeNumSameInNeighborhood + celltype_num_same;
        lineageNumSameInNeighborhood = lineageNumSameInNeighborhood + lineage_num_same;
        frameNumSameInNeighborhood = frameNumSameInNeighborhood + frame_num_same;
        
        %distance to nearest cell that is different from current (assuming it's in the neighborhood)
        %celltype
        for j = 1:length(inds)
            if embeddingStruct(inds(j)).celltype ~= embeddingStruct(i).celltype
                smallest_dist_to_other_celltype = smallest_dist_to_other_celltype + dists(j);
                break
            end
        end
        %lineage
        for j = 1:length(inds)
            if embeddingStruct(inds(j)).lineage ~= embeddingStruct(i).lineage
                smallest_dist_to_other_lineage = smallest_dist_to_other_lineage + dists(j);
                break
            end
        end
        %frame
        for j = 1:length(inds)
            if embeddingStruct(inds(j)).captureFrame ~= embeddingStruct(i).captureFrame
                smallest_dist_to_other_frame = smallest_dist_to_other_frame + dists(j);
                break
            end
        end
    end
end

%celltype 
%normalize to points tallied
celltypeNumSameInNeighborhood = celltypeNumSameInNeighborhood/num_tallied;
lineageNumSameInNeighborhood = lineageNumSameInNeighborhood/num_tallied;
frameNumSameInNeighborhood = frameNumSameInNeighborhood/num_tallied;

%add to output struct
scoreStruct3d.celltypeNumSameInNeighborhood = celltypeNumSameInNeighborhood;
scoreStruct3d.lineageNumSameInNeighborhood = lineageNumSameInNeighborhood;
scoreStruct3d.frameNumSameInNeighborhood = frameNumSameInNeighborhood;

%distance to nearest 'other'
%normalize by num tallied
smallest_dist_to_other_celltype = smallest_dist_to_other_celltype/num_tallied;
smallest_dist_to_other_lineage = smallest_dist_to_other_lineage/num_tallied;
smallest_dist_to_other_frame = smallest_dist_to_other_frame/num_tallied;

%normalize by average distance to any neighbor
smallest_dist_to_other_celltype = smallest_dist_to_other_celltype/dist_to_neighbor;
smallest_dist_to_other_lineage = smallest_dist_to_other_lineage/dist_to_neighbor;
smallest_dist_to_other_frame = smallest_dist_to_other_frame/dist_to_neighbor;

scoreStruct3d.celltypeSmallestDistToOther = smallest_dist_to_other_celltype;
scoreStruct3d.lineageSmallestDistToOther = smallest_dist_to_other_lineage;
scoreStruct3d.frameSmallestDistToOther = smallest_dist_to_other_frame;


%cluster compactness
%celltype
celltype_deviation = 0;
types_not_present = 1;
for j = 1:length(celltypes)
    celltype_current = celltypes(j);
    
    %grab indices,positions of cells of a given type
    celltype_indices = [];
    for i = 1:numPointsInEmbedding
        if embeddingStruct(i).celltype == celltype_current
            celltype_indices = [celltype_indices;i];
        end
    end
    celltype_points = points3d(celltype_indices,:);
    
    %prevent NaN values
    if isempty(celltype_points)
        types_not_present = types_not_present + 1;
        continue
    end
    
    %distances to centroid of celltype cluster
    centroid = [mean(celltype_points(:,1)),mean(celltype_points(:,2)),mean(celltype_points(:,3))];
    dists = pdist2(centroid,celltype_points);
    
    %log
    celltype_deviation = celltype_deviation + std(dists);
end
%normalize by number of organs
celltype_deviation = celltype_deviation/(length(celltypes)-types_not_present);


%normalize by same treatment on whole embedding
overall_centroid = [mean(points3d(:,1)),mean(points3d(:,2)),mean(points3d(:,3))];
overall_dists = pdist2(overall_centroid,points3d);
overall_std = std(overall_dists);
celltype_deviation = celltype_deviation/overall_std;

scoreStruct3d.celltypeCompactness = celltype_deviation;

%lineage
lineage_deviation = 0;
for j = 1:length(lineages)
    lineage_current = lineages(j);
    
    %grab indices,positions of cells of a given type
    lineage_indices = [];
    for i = 1:numPointsInEmbedding
        if embeddingStruct(i).lineage == lineage_current
            lineage_indices = [lineage_indices;i];
        end
    end
    lineage_points = points3d(lineage_indices,:);
    
    %distances to centroid of lineage cluster
    centroid = [mean(lineage_points(:,1)),mean(lineage_points(:,2)),mean(lineage_points(:,3))];
    dists = pdist2(centroid,lineage_points);
    
    %log
    lineage_deviation = lineage_deviation + std(dists);
end

%normalize by number of organs
lineage_deviation = lineage_deviation/length(lineages);

%normalize by same treatment on whole embedding
lineage_deviation = lineage_deviation/overall_std;

scoreStruct3d.lineageCompactness = lineage_deviation;


%frame
frame_deviation = 0;
for j = 1:length(distinct_frames)
    frame_current = distinct_frames(j);
    
    %grab indices,positions of cells of a given type
    frame_indices = [];
    for i = 1:numPointsInEmbedding
        if embeddingStruct(i).captureFrame == frame_current
            frame_indices = [frame_indices;i];
        end
    end
    frame_points = points3d(frame_indices,:);
    
    %distances to centroid of lineage cluster
    centroid = [mean(frame_points(:,1)),mean(frame_points(:,2)),mean(frame_points(:,3))];
    dists = pdist2(centroid,frame_points);
    
    %log
    frame_deviation = frame_deviation + std(dists);
end

%normalize by number of organs
frame_deviation = frame_deviation/length(distinct_frames);

%normalize by same treatment on whole embedding
frame_deviation = frame_deviation/overall_std;

scoreStruct3d.frameCompactness = frame_deviation;
end