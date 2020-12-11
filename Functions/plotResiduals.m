%% Function plotResidual
%%
% Inputs:
%   y_true  vector of true values
%   y_pred  vector of predicted/model values
%  optional
%    figure <n>  figure number
%    title <s>   title of figure
%    labels <s>
%    limits []   vector of lower and upper limits for signal
%%
% Output:
%   h   handle of created figure

% $Id$

function h = plotResiduals( y_true, y_pred, varargin )

p = inputParser;
p.addRequired('y_true',@isvector)
p.addRequired('y_pred',@isvector)
p.addParameter('figure',0,@isscalar)
p.addParameter('title','Residual plot',@ischar)
p.addParameter('labels',{'y_{true}','y_{pred}'},@iscell)
p.addParameter('limits',[],@ismatrix)
p.parse( y_true, y_pred, varargin{:} )
opts = p.Results;

if opts.figure > 0
    h = figure(opts.figure);
else
    h = figure;
end
clf

n = min( numel(y_true),numel(y_pred) );
y_true= y_true(1:n);
y_pred = y_pred(1:n);

subplot(4,5,[1:4,6:9,11:14,16:19])
plot( y_true, y_pred,'.')
xm=floor(min([y_true;y_pred]));
ym=ceil(max([y_true;y_pred]));
line( [xm,ym],[xm,ym],'Color','r','LineStyle','--')

if opts.limits, axis( opts.limits ), end
axis square

grid on
xlabel( opts.labels(1) )
ylabel( opts.labels(2) )

title( opts.title )

ec = ErrorCriteria( y_true, y_pred );
subplot(4,5,[5:5:20])
text( 0,0.9,sprintf( ['crit |  value   \n',...
    '-----+----------\n',...
    'MAE  | %7.3e\n',...
    'SSE  | %7.3e\n',...
    'MSE  | %7.3e\n',...
    'RMSE | %7.3e\n',...
    'NMSE | %7.3e\n',...
    'BFR  | %7.3e\n' ], ...
    ec.MAE,  ec.SSE,  ec.MSE, ...
    ec.RMSE, ec.NMSE, ec.BFR ), ...
    'FontName','Courier','VerticalAlignment','top' )
axis off
box off
title( 'Error criteria' )