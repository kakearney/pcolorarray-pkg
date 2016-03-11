function pcolorarray(varargin)
%PCOLORARRAY Display an array as a pcolor plot with value labels
%
% pcolorarray(var1, var2, ..., 'param', val, ...)
%
% This function displays an array as a table, with each cell colored
% according to the value of the cell.
%
% Input arguments:
%
%   var#:       two-dimensional double arrays to be displayed.  All
%               variables displayed will use a shared colormap. 
%
% Optional input arguments (passed as paramerer/value pairs, defaults in
% brackets) 
%
%   format:     formatting string used to label each value. ['%g']
%
%   cmap:       n x 3 colormap array [jet(64)]
%
%   nan:        string used to display NaNs. ['NaN']
%
%   cscale:     'linear' or 'log', specifies whether the color scaling for
%               the array uses a linear or logarithmic scale ['linear'] 
%
%   center0:    logical scalar, if true, center the colormap on 0 rather
%               than the mean of unique data [false]
%
%   rowname:    1 x nrow cell array of strings []
%
%   colname:    1 x ncol cell array of strings []
%
% You may also pass in any text properties except for 'units' and
% 'horizontalalignment'.  These properties will be applied to all displayed
% text.
	
% Copyright 2009 Kelly Kearney

%---------------------------
% Parse and check input
%---------------------------

% Parse input

Options = struct('format', '%g', ...
                 'cmap', jet(64), ...
                 'nan', 'NaN', ...
                 'cscale', 'linear', ...
                 'center0', false, ...
                 'rowname', [], ...
                 'colname', []);

isparam = cellfun(@ischar, varargin);
if any(isparam)
    idx = find(isparam,1);
    vars = varargin(1:idx-1);
    [Options, textprops] = parsepv(Options, varargin(idx:end), 'returnextra');
else
    vars = varargin;
    textprops = cell(0);
end


% Get names of input variables

nvar = length(vars);
for ivar = 1:nvar
	varnames{ivar} = inputname(ivar);
	if isempty(varnames{ivar})
		varnames{ivar} = sprintf('Input expression %d', ivar);
	end
end

% Check that variables are 2D and doubles 
% TODO allow other numeric formats if can be converted to double

for ivar = 1:nvar
	isgood = ndims(vars{ivar}) == 2 && strcmp('double', class(vars{ivar}));
end

if ~any(isgood)
	error('Can only display 2D doubles');
elseif any(~isgood)
	warning('Can only display 2D doubles; ignoring other input variables');
	vars = vars(isgood);
	nvar = length(vars);
end
    
% TODO: add checks for size if colname and rowname are supplied

hascollabel = ~isempty(Options.colname);
hasrowlabel = ~isempty(Options.rowname); 

if hasrowlabel
    Options.rowname = reshape(Options.rowname, [], 1);
end
if hascollabel
    Options.colname = reshape(Options.colname, 1, []);
end

%---------------------------
% Plot variables
%---------------------------

for iv = 1:nvar
    
    [nr,nc] = size(vars{iv});
   
    % Format for numbers
    
    vartext = arrayfun(@(x) sprintf(Options.format, x), vars{iv}, 'uni', 0);
    
    % Change display of NaNs

    [vartext{isnan(vars{iv})}] = deal(Options.nan);
    
    % Change row and column headers
    
    if hascollabel && ~hasrowlabel
        vars{iv} = [nan(1,nc); vars{iv}];
        vartext = [Options.colname; vartext];
        [nr,nc] = size(vars{iv});
    elseif hasrowlabel && ~hascollabel
        vars{iv} = [nan(nr,1) vars{iv}];
        vartext = [Options.rowname vartext]; 
        [nr,nc] = size(vars{iv});
    elseif hasrowlabel && hascollabel
        vars{iv} = [nan(1,nc+1); [nan(nr,1) vars{iv}]];
        vartext = [[{''} Options.colname]; [Options.rowname vartext]];
        [nr,nc] = size(vars{iv});
    end

    % Coordinates for each cell boundary and text label
    
    xcell = linspace(0,1,nc+1);
    ycell = linspace(0,1,nr+1);
    xtext = (xcell(1:end-1)+xcell(2:end))./2;
    ytext = (ycell(1:end-1)+ycell(2:end))./2;
    [xtext, ytext] = meshgrid(xtext, ytext);
    
    % Determine figure size needed to legibly show all cells

    hfigtemp = figure('visible', 'off');
    if verLessThan('matlab', '8.4.0')
        htxt = cellfun(@(x) text(0,0,x, 'units', 'pixels', textprops{:}), vartext);
        ext = cell2mat(get(htxt, 'extent'));
    else
        htxt = gobjects(size(vartext));
        for ii = 1:numel(htxt)
            htxt(ii) = text(0,0,vartext{ii},'units', 'pixels', textprops{:});
        end
        ext = cat(1, htxt.Extent);
    end
        
    delete(htxt);
    
    close(hfigtemp);

    cellsz = max(ext(:,3:4)) + 5;
    figsz = [cellsz(1)*nc cellsz(2)*nr];
    
	% Adjust values to log scale if necessary
	
	if strcmp(Options.cscale, 'log')
		vars{iv} = log10(vars{iv});
	end
	
	% Inf and -Inf ignored for color purposes
	
	vars{iv}(isinf(vars{iv})) = NaN;

	% Create figure
	
    hfig(iv) = figure('units', 'pixels', 'position', [100 100 figsz], 'name', varnames{iv}, 'menubar','none');
    hax(iv) = axes('position', [0 0 1 1]);
    pcolorpad(xcell, ycell, vars{iv});
    colormap(Options.cmap);
    
    cellfun(@(x,y,z) text(x,y,z, 'horiz', 'center', textprops{:}), num2cell(xtext), num2cell(ytext), vartext);
    
end

% Adjust all figures to use the same colormap limits

lims = cellfun(@minmax, vars, 'uni', 0);
lims = cat(1, lims{:});
lims = minmax(lims);
if Options.center0
	maxlim = max(abs(lims));
	lims = [-maxlim maxlim];
end
if diff(lims) == 0
    lims = lims + lims.*[-0.1 0.1];
end

set(hax, 'ydir', 'reverse', 'xtick', [], 'ytick', [], 'clim', lims);
