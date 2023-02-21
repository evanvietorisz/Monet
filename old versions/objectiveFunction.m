function [] = objectiveFunction(celltypes,organ_names,reducedfeature,reducedfeature3)
%objective function to score plots
%calculates how many of each cell's n nearest neighbors are the same type
%as it

%%%could still have errors. is also implemented badly. debug or delete

%number of neighbors to look at
n = 10;
num_organs = length(organ_names);

%first column: neighbor info for organ in 2d
%second column: number of cells of that organ
organs2d = zeros(num_organs,2);
%same in 3d
organs3d = zeros(num_organs,2);


for i = 1:length(reducedfeature)
    if celltypes(i) ~= 0
    
        [dists2d,inds2d] = pdist2(reducedfeature,reducedfeature(i,:),'euclidean','Smallest',n+1);
        inds2d = inds2d(2:end);
        
        [dists3d,inds3d] = pdist2(reducedfeature3,reducedfeature3(i,:),'euclidean','Smallest',n+1);
        inds3d = inds3d(2:end);

        for j = 1:n
            if celltypes(inds2d(j)) == celltypes(i)
                organs2d(celltypes(i),1) = organs2d(celltypes(i),1) + 1/n;
            end
        end
        organs2d(celltypes(i),2) = organs2d(celltypes(i),2) + 1;
        
        %reducedfeature3d is same length as reducedfeature
        for j = 1:n
            if celltypes(inds3d(j)) == celltypes(i)
                organs3d(celltypes(i),1) = organs3d(celltypes(i),1) + 1/n;
            end
        end
        organs3d(celltypes(i),2) = organs3d(celltypes(i),2) + 1;
    end
end

disp(organs2d)

for i = 1:num_organs
    organs2d(i,1) = organs2d(i,1)/organs2d(i,2);
    organs3d(i,1) = organs3d(i,1)/organs3d(i,2);
end

disp(organ_names)
disp('2d:')
disp(organs2d(:,1)')
disp('3d:')
disp(organs3d(:,1)')
end

