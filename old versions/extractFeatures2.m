function [featurevector] = extractFeatures2(embryo,starttime,endtime,k,sampling)
%given a loaded embryo/acetree data structure extract features for each
%timepoint concatenating them into one giant vector
%THIS VERSION ADJUSTED FOR MORE PARAMETERS

featurevector=[];
for t = starttime:sampling:endtime
    
    temp = computeFeaturesTimepoint2(embryo,k,t);
    featurevector = [featurevector;temp];
    
end
        
end

