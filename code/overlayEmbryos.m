function [embedding] = overlayEmbryos(title1,feat1,title2,feat2,neighborparam)
    %overlays two embryos in an embedding
    %   takes in featurevectors from two embryos (saved in the workspace),
    %   runs UMAP (2D) on both of them together, and plots the resultant
    %   embeddings colored by embryo
    
    %*requires umapFileEXchange to be in path
    
    %feats 1 and 2 are feature vectors from the workspace
    %titles 1 and 2 are names of the embryos being passed in (char array)
    %k is number of neighbors
    
    %compute embedding
    bigfeat = [feat1;feat2];
    embedding = run_umap(bigfeat,'n_neighbors',neighborparam);
    
    %color-code
    colors1 = repmat([1,0,0],size(feat1,1),1);
    colors2 = repmat([0,0,1],size(feat2,1),1);
    colors = [colors1; colors2];
    
    %plot
    figure('Name','Embedding','NumberTitle','off')
    scatter(embedding(:,1),embedding(:,2),10,colors)
    
    %decorate
    plottitle = ['Embedding of Two Embryos'];
    plottitle = [plottitle newline title1,' (red), ' title2, ' (blue)'];
    title(plottitle)
    xlabel('umap-1')
    ylabel('umap-2')
end
