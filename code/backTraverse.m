function [time_pt,ind_pt,ancestor_times,ancestor_inds] = backTraverse(start_time,start_ind,mode,time_parameter,embryo)
%traverses embryo dataset backwards a specified number of
%   generations/amount of time and returns times/indices of ancestors along
%   the way (times/indices correspond to last timestep at which the ancestor was itself)
%start_time is starting time
%start_ind is index of starting cell
%mode is 'generations' or 'timesteps', and specifies whether time_parameter
%refers to generations or timesteps
%embryo is embdat


time_pt = start_time;
ind_pt = start_ind;

counter = time_parameter;

ancestor_times = [];
ancestor_inds = [];

while counter > 0
    %walk back
    ind_pt = stepback(time_pt,ind_pt,embryo);
    time_pt = time_pt-1;
    %if crossed a cell division event, log information of progenitor,
    if embryo(time_pt).suc(ind_pt,2) ~= -1
        ancestor_times = [ancestor_times;time_pt];
        ancestor_inds = [ancestor_inds;ind_pt];
        
        %if generations, increment counter
        if mode == "generations"
            counter = counter-1;
        end
    end
    %if mode is timesteps, increment counter
    if mode == "timesteps"
        counter = counter-1;
    end
    
    %escape if at beginning of dataset
    if time_pt == 1
        break
    end
    
end
%if mode == generations, leaves you off at time/index of progenitor cell g generations ago
%at the last time it existed before dividing
%if mode==timesteps, leaves you off at that time
end

