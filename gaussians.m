function [val] = gaussians(x,y,z,sigma,derivative)
%computes value of specified distribution at point (x,y,z)
%   x,y,z are coordinates
%   sigma is stddev
%   derivative is what derivative you want
%   'zero' is 3d gaussian
%   'first_x', 'first_y', 'first_z' are NEGATIVE first derivatives wrt x,y,z directions
%   'laplacian' is NEGATIVE laplacian


if derivative == "zero"
    val = 1/((2*pi)^(3/2)*sigma^3)*exp(-(x.^2+y.^2+z.^2)/(2*sigma^2));
end

if derivative == "first_x"
    val = x/(sigma^2).*(1/((2*pi)^(3/2)*sigma^3)*exp(-(x.^2+y.^2+z.^2)/(2*sigma^2)));
end

if derivative == "first_y"
    val = y/(sigma^2).*(1/((2*pi)^(3/2)*sigma^3)*exp(-(x.^2+y.^2+z.^2)/(2*sigma^2)));
end

if derivative == "first_z"
    val = z/(sigma^2).*(1/((2*pi)^(3/2)*sigma^3)*exp(-(x.^2+y.^2+z.^2)/(2*sigma^2)));
end

if derivative == "laplacian"
    val = (-(x.^2+y.^2+z.^2)/sigma^4-3/sigma^2).*(1/((2*pi)^(3/2)*sigma^3)*exp(-(x.^2+y.^2+z.^2)/(2*sigma^2)));
end

end

