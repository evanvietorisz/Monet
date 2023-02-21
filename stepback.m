function prev_ind = stepback(time,ind,embryo)
%returns index of cell in prev time that tracks to cell 'ind' at current time

%embryo is embdat_stabilized
%time, ind are time and index of cell in embdat(time).finalpoints you want
%to start from

%note: behavior/behavior of scripts that call this may be not throw and error
%for datasets with tracking errors. be sure to double check tracking before
%using this function and backTraverse/forwardTraverse. 


if time == 1
    prev_ind = ind;
    return
end

found_it = false;

for count = 1:size(embryo(time-1).suc,1)
    if embryo(time-1).suc(count,1) == ind || embryo(time-1).suc(count,2) == ind
        prev_ind = count;
        found_it = true;
        break
    end
end

end