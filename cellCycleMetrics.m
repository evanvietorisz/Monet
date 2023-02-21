function [beginning_of_current_cycle,end_of_current_cycle,duration_of_current_cell_cycle,time_until_next_division,...
    time_since_last_division,fraction_of_current_cycle_elapsed,time_of_next_division] = cellCycleMetrics(t,i,embryo)
%returns metrics of cell cycle for cellin question -- embdat(t).finalpoints(i,:)
%   beginning_of_current_cell_cycle (absolute time)
%   end_of_cell_current_cycle (absolute time)
%   duration_of_current_cell_cycle
%   time_until_next_division
%   time_since_last_division
%   fraction_of_current_cycle_elapsed


%debugging
%['time = ',num2str(t),', i = ',num2str(i)]
%t, i are time and index of cell you want to start from
%embryo is embdat

%extract beginning and end times of current cell cycle

%beginning (inclusive)
ind = i;
time = t;

if time ~= 1
    
    while embryo(time-1).suc(stepback(time,ind,embryo),2) == -1
        
        
        ind = stepback(time,ind,embryo);
        time = time-1;
        if time == 1
            break
        end
    end
end

beginning_of_current_cycle = time;

%end (inclusive)
ind = i;
time = t;
endtime = length(embryo)-1;
if time ~= endtime
    while embryo(time).suc(ind,2) == -1
        ind = embryo(time).suc(ind,1);
        time = time + 1;

        %escape if reach end of embryo
        if time == endtime
            break
        end
        
        %escape if cell dies, treat as site of its next division
        if ind == -1
            break
        end
        
        
    end
end
end_of_current_cycle = time;


duration_of_current_cell_cycle = end_of_current_cycle - beginning_of_current_cycle + 1;

time_until_next_division = end_of_current_cycle - t;
time_since_last_division = t - beginning_of_current_cycle + 1;
fraction_of_current_cycle_elapsed = (t - beginning_of_current_cycle + 1)/(end_of_current_cycle - beginning_of_current_cycle + 1);
time_of_next_division = end_of_current_cycle;%antiquated

%note
%beginning_of_current_cell_cycle is first time at which that cell is what
%it is (is no longer the cell that came before it)
%end_of_cell_current_cycle is the last time at which it still is, not the
%time of the next division (= end_of_cell_current_cycle+1)
%this is the source of the +1 convention in defining the other quantities













end

