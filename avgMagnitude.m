function [normed_array] = avgMagnitude(array)
%takes in an array where each row is a vector or multiple vectors,
%returns the array normalized by the average magnitude of the vectors

column_of_vecs = zeros(size(array,1)*size(array,2)/3,3);

c=1;
for i = 1:size(array,1)
    for j = 1:3:size(array,2)
        column_of_vecs(c,:) = array(i,j:j+2);
        c = c+1;
    end
end

column_of_magnitudes = zeros(size(column_of_vecs,1),1);
for i = 1:size(column_of_vecs,1)
    column_of_magnitudes(i) = norm(column_of_vecs(i,:));
end

avg_mag = mean(column_of_magnitudes);

normed_array = array/avg_mag;

end

