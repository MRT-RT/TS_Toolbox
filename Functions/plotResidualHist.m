%% Function plotResidualHist
%%
% Inputs:
%   y_true  vector of true values
%   y_pred  vector of predicted/model values
%  optional:
%    nbins <n>      number of histogram classes
%    figure <n>     figure number
%    title <s>      title of figure
%    labels <s>     label x/y axis 
%    limits []      vector of lower and upper limits for signal
% Output:
%   h   handle of created figure

% $Id$

function h = plotResidualHist( y_true, y_pred, varargin )

p = inputParser;
p.addRequired('y_true',@isvector)
p.addRequired('y_pred',@isvector)
p.addParameter('nbins',10,@isscalar)
p.addParameter('figure',0,@isscalar)
p.addParameter('title','Residual histogram',@ischar)
p.addParameter('labels',{'residuals','frequency'},@iscell)
p.addParameter('limits',[],@ismatrix)
p.parse( y_true, y_pred, varargin{:} )
opts = p.Results;

if opts.figure > 0
    h = figure( opts.figure );
else
    h = figure;
end
clf

% Residual
r = y_true - y_pred;
nr = numel( r );
qr = myQuantile( r );
mr = mean( r );
sr = std( r );

histogram( r, opts.nbins )
grid on

set(gca ,'FontSize',14 )
xlabel( opts.labels(1),'FontSize',14  )
ylabel( opts.labels(2),'FontSize',14  )
title( opts.title,'FontSize',14 )

% Limits for histogramm (only x or x&y)
if opts.limits
    l = length( opts.limits );
    if l == 2
        xlim( opts.limits )
    elseif l == 4
        axis( opts.limits )
    end
end


s = sprintf( [ 'Statistiscs:\n\n', ...
    'N = %d data points\n\n', ...
    'crit   |  value   \n',...
    '-------+-----------\n',...
    'mean   | %+7.4e\n',...
    'std    | %+7.4e\n',...
    'min    | %+7.4e\n',...
    'q25    | %+7.4e\n',...
    'median | %+7.4e\n',...
    'q75    | %+7.4e\n',...
    'max    | %+7.4e\n'],nr,mr,sr,qr);
xl = get( gca, 'Xlim' );
dx = diff( xl );
yl = get( gca, 'Ylim' );
dy = diff( yl );

text( xl(2)-dx/5,yl(2)-dy/6, s, ...
    'FontName','Courier','FontSize',10,'VerticalAlignment','top' )

h.PaperOrientation = 'landscape';
h.WindowState = 'maximized';
