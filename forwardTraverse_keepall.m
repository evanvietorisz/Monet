function [cells_to_look_at,time_pt,all_later_cell_times,all_later_cell_inds] = forwardTraverse_keepall(start_time,cells_to_look_at,mode,time_parameter,embryo)
%same as forwardTraverse but keeps a log of all cells traversed over


%traverse forward through embryo, keeping track of the descendants of a
%   predetermined list of cells 
%start_time is time to start from 
%cells_to_look_at is indices of cells at start_time whose descendants you
%want to trace
%embryo is embdat

emb_size = length(embryo);

time_pt = start_time;



%indices between these arrays correspond to each other
all_later_cell_times = [];
all_later_cell_inds = [];




if mode == "timesteps" 
    while time_pt < time_parameter
        cells_to_look_at_next = [];
        for j = 1:length(cells_to_look_at)

            descendants = embryo(time_pt).suc(cells_to_look_at(j),:);
            if descendants(2) == -1
                cells_to_look_at_next = [cells_to_look_at_next,descendants(1)];
            else
                cells_to_look_at_next = [cells_to_look_at_next,descendants];
            end
        end
        cells_to_look_at = cells_to_look_at_next;
        %remove any dead cells that got in
        cells_to_look_at(find(cells_to_look_at==-1)) = [];
        
        time_pt = time_pt+1;
        
        %add to cells visited
        all_later_cell_inds = [all_later_cell_inds;cells_to_look_at'];
        all_later_cell_times = [all_later_cell_times;repmat([time_pt],size(cells_to_look_at',1),1)];
        

        if time_pt >= emb_size-1
            break
        end
        
        %no need for escape since user wouldn't put in an end time beyond
        %the duration of the dataset
    end
end

if mode == "generations"
    generations_passed = 0;
    while generations_passed < time_parameter
        cells_to_look_at_next = [];
        
        %if cell dies, break
        if embryo(time_pt).suc(cells_to_look_at(1)) == -1
            break
        end
   
        for j = 1:length(cells_to_look_at)

            descendants = embryo(time_pt).suc(cells_to_look_at(j),:);
            if descendants(2) == -1
                cells_to_look_at_next = [cells_to_look_at_next,descendants(1)];
            else
                cells_to_look_at_next = [cells_to_look_at_next,descendants];
            end
        end
        cells_to_look_at = cells_to_look_at_next;
        time_pt = time_pt+1;
 
        %add to cells visited
        all_later_cell_inds = [all_later_cell_inds;cells_to_look_at'];
        all_later_cell_times = [all_later_cell_times;repmat([time_pt],size(cells_to_look_at',1),1)];
        
        
        
        %log generations passed
        if descendants(2) ~= -1
            generations_passed = generations_passed + 1;
        end
        
        %escape if at end of embryo
        if time_pt >= emb_size-1
            break
        end
    end
    
end


%leaves you off with time_pt = end_time and cells_to_look_at = the indices
%of the decscendants of the first cells_to_look_at at time end_time



end