function [frac] = fracPointingOut(points,normals)
%fracPointingOut Returns fraction of points whose calculated normals point
%   away from the origin

n = size(points,1);

%figure out how normals oriented
points_plus_norms = points + normals;
orientations = zeros(num_points,1);
for i = 1:n
    %if normal points inward, set value to 1
    if norm(points(i)) >= norm(points_plus_norms(i))
        orientations(i) = 1;
    end
end

%fraction of normals oriented out
frac = 0;
for i = 1:n
    frac = frac + orientations(i);
end
frac = frac/num_points;

end

