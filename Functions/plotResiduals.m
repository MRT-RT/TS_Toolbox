%% Function plotResiduals
%%
% Inputs:
%   y_obsv  vector of observed values
%   y_pred  vector of predicted/model values
%  optional:
%    parameter <n>  number of mdel parameters
%    figure <n>     figure number
%    title <s>      title of figure
%    labels <s>     label x/y axis 
%    limits []      vector of lower and upper limits for signal
% Output:
%   h   handle of created figure

% $Id$

function h = plotResiduals( y_obsv, y_pred, varargin )

p = inputParser;
p.addRequired('y_obsv',@isvector)
p.addRequired('y_pred',@isvector)
p.addParameter('parameter',0,@isscalar)
p.addParameter('figure',0,@isscalar)
p.addParameter('title','Residual plot',@ischar)
p.addParameter('labels',{'y_{obsv}','y_{pred}'},@iscell)
p.addParameter('limits',[],@ismatrix)
p.parse( y_obsv, y_pred, varargin{:} )
opts = p.Results;

if opts.figure > 0
    h = figure(opts.figure);
else
    h = figure;
end
clf

n = min( numel(y_obsv),numel(y_pred) );
y_obsv= y_obsv(1:n);
y_pred = y_pred(1:n);

%subplot(4,5,[1:4,6:9,11:14,16:19])
plot( y_obsv, y_pred,'.','MarkerSize',10)
xm=floor(min([y_obsv;y_pred]));
ym=ceil(max([y_obsv;y_pred]));
line( [xm,ym],[xm,ym],'Color','r','LineStyle','--','LineWidth',2)

if opts.limits, axis( opts.limits ), end
axis square

grid on
set(gca ,'FontSize',14 )
xlabel( opts.labels(1),'FontSize',14  )
ylabel( opts.labels(2),'FontSize',14  )

title( opts.title,'FontSize',14 )

ec = ErrorCriteria( y_obsv, y_pred, opts.parameter );

s = sprintf( [ '\n Error criteria:\n\n', ...
    ' N = %d data points\n\n', ...
    ' crit |  value   \n',...
    ' -----+-----------\n',...
    ' MAE  | %7.4e\n',...
    ' SSE  | %7.4e\n',...
    ' MSE  | %7.4e\n',...
    ' RMSE | %7.4e\n',...
    ' NMSE | %7.4e\n',...
    ' BFR  | %7.4e\n' ], ...
    n,ec.MAE,  ec.SSE,  ec.MSE, ...
    ec.RMSE, ec.NMSE, ec.BFR );
if opts.parameter > 0
    s = [ s, sprintf( [ ...
        ' AIC  | %7.4e\n',...
        ' BIC  | %7.4e\n',...
        ' nPar | %4d\n' ], ...
        ec.AIC, ec.BIC,opts.parameter ) ];
end
xl=get(gca,'Xlim');
dx=diff(xl);
yl=get(gca,'Ylim');

%text( xl(2)-dx/4,yl(1), s, ...
%    'FontName','Courier','FontSize',10,'VerticalAlignment','bottom' )
text( xl(2),yl(1), s, ...
    'FontName','Courier','FontSize',10,'VerticalAlignment','bottom',...
    'Units','normalized','Position',[0.7,0.05])

