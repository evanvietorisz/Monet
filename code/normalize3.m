function [P] = normalize3(percellfeature,feature_names_list,sizes_list,norm_types_list,embryo)
%function that takes in featurevector and normalizes each feature according
%to prescribed method

%make copy of featurevector P = normalized percellfeature
P = percellfeature;


%get indices of each feature from sizes
indices_list = zeros(size(feature_names_list,1),2);
c = 1;
for i = 1:length(feature_names_list)
    indices_list(i,1) = c;
    indices_list(i,2) = c + sizes_list(i)-1;
    c = c + sizes_list(i);
end


%for later

%obtain first time at which emb is 4 cells
fourcelltime = 1;
while size(embryo(fourcelltime).finalpoints,1) < 4
    fourcelltime = fourcelltime + 1;
end

%obtain last time at which there is an entry for timedata (movement starts)
lasttime = length(embryo) - 1; %just a consequence of how the emb is loaded, assuming it is loaded up until movement starts

%obtain max number of cells emb has before movement starts
maxnumcells = size(embryo(lasttime).finalpoints,1);


%normalize each feature
for i = 1:size(feature_names_list,1)
    

    if norm_types_list(i) == "scale_4cell_to_end"
        %linearly scale all entries on [0,1] where 0 is when the emb is 4
        %cells and 1 is the end of the dataset (when the worm starts moving)
        P(:,indices_list(i,1):indices_list(i,2)) = rescale(P(:,indices_list(i,1):indices_list(i,2)),0,1,'InputMin',fourcelltime,'InputMax',lasttime);

    elseif norm_types_list(i) == "scale_0_to_max_num_cells"
        %linearly scale all entries on [0,1] where 0 is 4 cells (minimum number of cells if you begin counting time at the 4 cell stage) and
        %1 is the number of cells at end of the dataset (when the worm starts moving)
        P(:,indices_list(i,1):indices_list(i,2)) = rescale(P(:,indices_list(i,1):indices_list(i,2)),0,1,'InputMin',4,'InputMax',maxnumcells);
        
    elseif norm_types_list(i) == "column_zscore"
        %take zscore of column ( indices_list(i,1)==indices_list(i,2) )
        P(:,indices_list(i,1):indices_list(i,2)) = zscore(P(:,indices_list(i,1):indices_list(i,2)));
        
    elseif norm_types_list(i) == "all_zscore"
        %take zscore of all elements
        P(:,indices_list(i,1):indices_list(i,2)) = zscore(P(:,indices_list(i,1):indices_list(i,2)),0,'all');
        
    elseif norm_types_list(i) == "divide_avg_mag"
        %divide array elements by avg magitude of a vector in the array
        P(:,indices_list(i,1):indices_list(i,2)) = avgMagnitude(P(:,indices_list(i,1):indices_list(i,2)));
    end
    
end

end
