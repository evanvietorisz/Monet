function [featurevector] = extractFeatures3(embryo,starttime,endtime,sampling,feats_and_norms,param_values)
%given a loaded embryo/acetree data structure extract features for each
%timepoint concatenating them into one giant vector

featurevector=[];
for t = starttime:sampling:endtime
    
    temp = computeFeaturesTimepoint3(embryo,t,feats_and_norms,param_values);
    featurevector = [featurevector;temp];
    
end
        
end
