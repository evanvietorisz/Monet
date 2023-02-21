function [valid_indices] = embeddingLookup(embeddingStruct,xs,ys)
%convenient way to look at what is going on in the embeddings
%note: only for 2d since that is more useful
%embeddingStruct is the embeddingStruct of interest
%xs is array containing [xmin,xmax]
%ys is array containing [ymin,ymax]

valid_indices = [];

c = 1;
for i = 1:size(embeddingStruct,2)
    if xs(1) < (embeddingStruct(i).twoDpoint(1)) && (embeddingStruct(i).twoDpoint(1) < xs(2)) &&...
            (ys(1) < embeddingStruct(i).twoDpoint(2)) && (embeddingStruct(i).twoDpoint(2) < ys(2))
        valid_indices(c) = i;
        c = c + 1;
    end
end

end