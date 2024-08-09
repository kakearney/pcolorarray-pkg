function h = pcolorpad(varargin)
%PCOLORPAD Pseudocolor plot defined by edges of cells
%
% h = pcolorpad(c)
% h = pcolorpad(x,y,c)
% h = pcolorpad(ax, ...)
% 
% This function is nearly identical to pcolor, except it doesn't drop the
% last row and column of an array.
%
% Input variables:
%
%   c:  m x n array of color values to be plotted as a surface
%
%   x:  x-coordinate data for the surface, defining the x-coordinates of
%       the surface cell vertices.  Can be either an 1 x n vector, 1 x n+1
%       vector, m x n array, m x n+1 array, m+1 x n array, or m+1 x n+1
%       array.  If the number of columns in x is one more than the number
%       of columns in c, then the last column of the c array will be
%       visible.  If x is a 2D array, its size must match that of y.
%
%   y:  y-coordinate data for the surface, defining the y-coordinates of
%       the surface cell vertices.  Can be either an m x 1 vector, m+1 x 1
%       vector, m x n array, m x n+1 array, m+1 x n array, or m+1 x n+1
%       array.  If the number of rows in y is one more than the number
%       of rows in c, then the last row of the c array will be
%       visible.  If y is a 2D array, its size must match that of x.
%  
%   ax: handle to axis where pcolor surface is plotted
%
% Optional input variables:
%
%   all parameter/value pairs available for a pcolor or surf plot can be
%   passed to this function.  
%
% Output variables:
%
%   h:  handle to plotted surface object

% Copyright 2017 Kelly Kearney

% Look for axis input

isax = cellfun(@(x) isscalar(x) && ishandle(x) && strcmp(get(x, 'type'), 'axes'), varargin);
if any(isax)
    ax = varargin{isax};
    args = varargin(~isax);
else
    args = varargin;
    ax = gca;
end

% Look for surf/pcolor parameter-value pairs

isstr = cellfun(@ischar, args);
if any(isstr)
    idx = find(isstr,1);
    params = args(idx:end);
    xyz = args(1:idx-1);
else
    xyz = args;
    params = {};
end

% Check input

if length(xyz) == 1
    z = xyz{1};
    validateattributes(z, {'numeric'}, {'2d'}, 'pcolorpad', 'c');
    
    x = 1:size(z,2)+1;
    y = 1:size(z,1)+1;
    
elseif length(xyz) == 3
    x = xyz{1};
    y = xyz{2};
    z = xyz{3};
    
    validateattributes(x, {'numeric','datetime'}, {'2d'}, 'pcolorpad', 'x'); % TODO: make as flexible as pcolor
    validateattributes(y, {'numeric'}, {'2d'}, 'pcolorpad', 'y');
    validateattributes(z, {'numeric'}, {'2d'}, 'pcolorpad', 'c');
    
else
    error('Incorrect number of input arguments');
end

% Check sizes of x and y to determine where to pad

[ny, nx] = size(z);

if isvector(x)
    x = reshape(x, 1, []);
end
if isvector(y)
    y = y(:);
end

if ~ismember(size(x,2), [nx nx+1])
    error('Number of columns in x array must match size(c,2) or size(c,2)+1');
end
if ~ismember(size(y,1), [ny ny+1])
    error('Number of rows in y array must match size(c,1) or size(c,1)+1');
end

% Pad x

if size(x,2) == nx+1
    z = [z nan(ny,1)];
end 

% Pad y

if size(y,1) == ny+1
    z = [z; nan(1, size(z,2))];
end

% Plot
    
h = pcolor(ax, x, y, z, params{:});
