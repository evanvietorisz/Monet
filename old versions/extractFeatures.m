function [ featurevector ] = extractFeatures( embryo,starttime,endtime,k,sampling )
%given a loaded embryo/acetree data structure extract features for each
%timepoint concatenating them into one giant vector


featurevector=[];
for i=starttime:sampling:endtime
    temp=computeFeaturesTimepoint(embryo(i).finalpoints,k);
featurevector=[featurevector;temp];
end
        


end

