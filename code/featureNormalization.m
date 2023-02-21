function [treatedFeatureVector] = featureNormalization(featurevector,k,embdat_stabilized,starttime,sampling,endtime)
%preprocessing for large feature vector
%normalizes each column of dataset separately
%makes mean = 0, std = 1, vals ~~~ [-3,3]

%manual for large feature vector

numCols = size(featurevector,2);


%{
%all except rel position
componentwise = [1,k;4*k+1,5*k;8*k+1,9*k;9*k+2,10*k+1];
scalars = [9*k+1];
vectorquantities = [k+1,4*k;5*k+1,8*k];
%{
%all
componentwise = [1,k;4*k+1,5*k;8*k+1,9*k;9*k+5,10*k+4];
scalars = [9*k+4];
vectorquantities = [k+1,4*k;5*k+1,8*k;9*k+1,9*k+3];
%}
%}

%local only
componentwise = [1,k;4*k+1,5*k;8*k+1,9*k;9*k+1,10*k];
scalars = [];
vectorquantities = [k+1,4*k;5*k+1,8*k];

%{
%looking for branching pattern

%distfeat', rel_position, num_cells
componentwise = [1,k];
scalars = [k+4];
vectorquantities = [k+1,k+3];
%}

%{
%whole time point divided through by average inter-cell distance
for i = 1:size(featurevector,1)
    for j = 1:size(componentwise,1)
        %whole time point divided through by average inter-cell distance
        %distfeat is 1:kin each row
        featurevector(i,:) = (featurevector(i,:))/mean(featurevector(i,1:k));
    end
end
%}

%{
%normalizing the local features/everything that depends on them by their
%mean value
for i = 1:size(featurevector,1)
    for j = 1:size(componentwise,1)
        %each feature individually
        featurevector(i,componentwise(j,1):componentwise(j,2)) = (featurevector(i,componentwise(j,1):componentwise(j,2)))/(mean(featurevector(i,componentwise(j,1):componentwise(j,2))));
    end
end
%}

%{
%each vector feature divided through by the average magnitude of each
%vector
for i = 1:size(featurevector,1)
    for j = 1:size(vectorquantities,1)
        horizontal = featurevector(i,vectorquantities(j,1):vectorquantities(j,2));
        columnvec = [];
        for h = 1:3:length(horizontal)
            columnvec = [columnvec;horizontal(h:h+2)];
        end
        avgMagnitude = norm(mean(columnvec));
        normedColumnvec = columnvec/avgMagnitude;
        %elongate normedColumnvec
        normedHorizontal = [];
        for h = 1:size(normedColumnvec)
            normedHorizontal = [normedHorizontal,normedColumnvec(h,:)];
        end
        featurevector(i,vectorquantities(j,1):vectorquantities(j,2)) = normedHorizontal;
    end
end
%}

%{
%zscoring each feature separately
for i = 1:size(featurevector,1)
    for j = 1:size(componentwise,1)
        featurevector(i,componentwise(j,1):componentwise(j,2)) = zscore(featurevector(i,componentwise(j,1):componentwise(j,2)));
    end
end


for i = 1:size(featurevector,1)
    for j = 1:size(vectorquantities,1)
        horizontal = featurevector(i,vectorquantities(j,1):vectorquantities(j,2));
        columnvec = [];
        for h = 1:3:length(horizontal)
            columnvec = [columnvec;horizontal(h:h+2)];
        end
        avgvec = mean(columnvec);
        avgMagnitude = norm(avgvec);

        meanCenteredColumnvec = columnvec - avgvec;
        
        stdMagnitude = std(vecnorm(meanCenteredColumnvec,2,2));
        
        zscoredColumnvec = meanCenteredColumnvec/stdMagnitude;
        
        %elongate normedColumnvec
        normedHorizontal = [];
        for h = 1:size(zscoredColumnvec)
            normedHorizontal = [normedHorizontal,zscoredColumnvec(h,:)];
        end
        featurevector(i,vectorquantities(j,1):vectorquantities(j,2)) = normedHorizontal;
    end
end
%}

%treatedFeatureVector = featurevector;
%return

%{
%dividing through each feature separately
for j = 1:size(componentwise,1)
    for i = componentwise(j,1):componentwise(j,2)
        featurevector(:,i) = (featurevector(:,i)-mean(featurevector(:,i)))/std(featurevector(:,i));
    end
end
%}

%timepoint normalizations for anthony 10/13/20
%{
%time-block normalization by average inter-cell distance
c = 1;
for t = starttime:sampling:endtime
    n = length(embdat_stabilized(t).finalpoints);
    %divide through by distance feature (elements 1:k)
    featurevector(c:c+n-1,:) = featurevector(c:c+n-1,:)/mean(featurevector(c:c+n-1,1:k),'all');
    c = c + n;
end
%}


%{
%time-block normalization by average value (magnitude for vector features) for each feature
c = 1;
for t = starttime:sampling:endtime
    n = length(embdat_stabilized(t).finalpoints);
    %divide through by feature average/magnitude
    for i = 1:size(componentwise,1)
        featurevector(c:c+n-1,componentwise(i,1):componentwise(i,2)) = featurevector(c:c+n-1,componentwise(i,1):componentwise(i,2))/mean(featurevector(c:c+n-1,componentwise(i,1):componentwise(i,2)),'all');
    end
    c = c + n;
end

c = 1;
for t = starttime:sampling:endtime
    n = length(embdat_stabilized(t).finalpoints);
    
    for i = 1:size(vectorquantities,1)
        avg_mags = [];
        for j = c:c+n-1
            horizontal = featurevector(j,vectorquantities(i,1):vectorquantities(i,2));
            column = [];
            for l = 1:3:length(horizontal)
                column = [column;horizontal(l:l+2)];
            end
            avg_mag = 0;
            for m = 1:size(column,1)
                avg_mag = avg_mag + norm(column(m,:))/size(column,1);
            end
            avg_mags = [avg_mags;avg_mag];
        end
        average_magnitude = mean(avg_mags);
        featurevector(c:c+n-1,vectorquantities(i,1):vectorquantities(i,2)) = featurevector(c:c+n-1,vectorquantities(i,1):vectorquantities(i,2))/average_magnitude;
    end
    c = c + n;
end
%}

%{
%1/21/20
%distfeat (per timeblock, zscore on whole thing), displacement (per timeblock meancenter,divide by std
%of all values), developmental time (NOT per timeblock, zscore of whole column)

componentwise = [1,k];
scalars = [k+4];
vectorquantities = [k+1,k+3];

c = 1;
for t = starttime:sampling:endtime
    n = length(embdat_stabilized(t).finalpoints);
    %distfeat
    featurevector(c:c+n-1,1:k) = zscore(featurevector(c:c+n-1,1:k),0,'all');
    
    %displacementvec
    featurevector(c:c+n-1,k+1:k+3) = (featurevector(c:c+n-1,k+1:k+3) - mean(featurevector(c:c+n-1,k+1:k+3),1))/std(featurevector(c:c+n-1,k+1:k+3),0,'all');
    
    c = c + n;
end

%time
featurevector(:,k+4) = zscore(featurevector(:,k+4));
%}

%{
%1/21/20
%distfeat (per timeblock, zscore on whole thing),developmental time (NOT per timeblock, zscore of whole column)
componentwise = [1,k];
scalars = [k+1];
c = 1;
for t = starttime:sampling:endtime
    n = length(embdat_stabilized(t).finalpoints);
    %distfeat
    featurevector(c:c+n-1,1:k) = zscore(featurevector(c:c+n-1,1:k),0,'all');
    c = c + n;
end
featurevector(:,k+1) = zscore(featurevector(:,k+1));
%}

%{
%10/21/20
%distfeat (per timeblock, zscore on whole thing), displacement (per timeblock meancenter,divide by std
%of all values)

componentwise = [1,k];
vectorquantities = [k+1,k+3];

c = 1;
for t = starttime:sampling:endtime
    n = length(embdat_stabilized(t).finalpoints);
    %distfeat
    featurevector(c:c+n-1,1:k) = zscore(featurevector(c:c+n-1,1:k),0,'all');
    
    %displacementvec
    featurevector(c:c+n-1,k+1:k+3) = (featurevector(c:c+n-1,k+1:k+3) - mean(featurevector(c:c+n-1,k+1:k+3),1))/std(featurevector(c:c+n-1,k+1:k+3),0,'all');
    
    c = c + n;
end
%}
 
%{
%10/21/20
%distfeat (per timeblock, zscore on whole thing), displacement (per timeblock meancenter,divide by std
%of all values)

componentwise = [1,k];

c = 1;
for t = starttime:sampling:endtime
    n = length(embdat_stabilized(t).finalpoints);
    %distfeat
    featurevector(c:c+n-1,1:k) = zscore(featurevector(c:c+n-1,1:k),0,'all');
    c = c + n;
end
%}

%{
%10/21/20
%distfeat (per timeblock, zscore on whole thing), displacement (per timeblock meancenter,divide by std
%of all values)

vectorquantities = [1,3];

c = 1;
for t = starttime:sampling:endtime
    n = length(embdat_stabilized(t).finalpoints);
    %displacementvec
    featurevector(c:c+n-1,1:3) = (featurevector(c:c+n-1,1:3) - mean(featurevector(c:c+n-1,1:3),1))/std(featurevector(c:c+n-1,1:3),0,'all');

    c = c + n;
end
%}

%10/21/20
%distfeat (per timeblock, zscore on whole thing), displacement (per timeblock meancenter,divide by std
%of all values)

vectorquantities = [1,3];
scalars = [4]

c = 1;
for t = starttime:sampling:endtime
    n = length(embdat_stabilized(t).finalpoints);
    %displacementvec
    featurevector(c:c+n-1,1:3) = (featurevector(c:c+n-1,1:3) - mean(featurevector(c:c+n-1,1:3),1))/std(featurevector(c:c+n-1,1:3),0,'all');

    c = c + n;
end
featurevector(:,4) = zscore(featurevector(:,4));





treatedFeatureVector = featurevector;








return


%{
%normal version
%z-score
for j = 1:size(componentwise,1)
    for i = componentwise(j,1):componentwise(j,2)
        featurevector(:,i) = (featurevector(:,i)-mean(featurevector(:,i)))/std(featurevector(:,i));
    end
end


%vector ones that need three component treatment
for j = 1:size(vectorquantities,1)
    featline = [];
    for i = vectorquantities(j,1):vectorquantities(j,2)
        featline = [featline;featurevector(i,:)'];
    end
    dev = std(featline);
    for i = vectorquantities(j,1):vectorquantities(j,2)
        featurevector(:,i) = (featurevector(:,i)-mean(featurevector(:,i)))/std(featline);
    end
end

if ~isempty(scalars)
    for j = 1:length(scalars)
        i = scalars(j);%note no bounds; just a single value
        featurevector(:,i) = (featurevector(:,i)-mean(featurevector(:,i)))/std(featurevector(:,i));
    end
end


treatedFeatureVector = featurevector;
end
%}
