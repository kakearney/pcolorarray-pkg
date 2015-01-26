function h = pcolorpad(x,y,z, varargin)
%PCOLORPAD Pseudocolor plot defined by edges of cells
%
% h = pcolorpad(x,y,z)
% 
% This function is nearly identical to pcolor, except it doesn't drop the
% last row and column of an array.
%
% Input variables:
%
%   x:  n x 1 array defining edges in x direction
%
%   y:  m x 1 array defining edges in y direction
%
%   z:  m-1 x n-1 array of z values

if nargin == 1
    z = x;
    x = 1:(size(z,2)+1);
    y = 1:(size(z,1)+1);
end

if isequal(size(z), [length(y) length(x)])
    x2 = [x(:); NaN];
    y2 = [y(:); NaN];
elseif isequal(size(z), [length(y)-1 length(x)-1])
    x2 = x;
    y2 = y;
elseif isequal(size(x), size(y), size(z)+[1 1])
    x2 = x;
    y2 = y;
else
    error('z must be length(y) x length(x) or length(y)-1 x length(x)-1 array');
end

[nrow, ncol] = size(z);
z2 = nan(nrow+1, ncol+1);
z2(1:nrow, 1:ncol) = z;

h = pcolor(x2, y2, z2, varargin{:});
