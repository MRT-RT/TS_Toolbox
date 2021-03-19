%% Function plotRuleActivation
%%
% Inputs:
%   u      input vector/maxtrix
%   y      output vector
%   model  TS Model or msf(u,y)
%  optional:
%    parameter <n>  number of mdel parameters
%    figure <n>     figure number
%    title <s>      title of figure
%    labels <s>     label t/u_i/y axis
%    limits []      vector of lower and upper limits for signals
% Output:
%   h   handle of created figure

% $Id$

function h = plotRuleActivation( u, y, model, varargin )


valModel = @(x) isa(x,'TSModel') || ismatrix(x);

p = inputParser;
p.addRequired('u',@ismatrix)
p.addRequired('y',@isvector)
p.addRequired('model',valModel)
p.addParameter('figure',0,@isscalar)
p.addParameter('labels',{'time','input u','output y'},@iscell)
p.addParameter('limits',[],@ismatrix)
p.addParameter('title','Rule activation',@ischar)
p.parse( u,y, model, varargin{:} )
opts = p.Results;

% TS model or mu given?
if isa( opts.model, 'TSModel' )
    mu = getMSF( model );
elseif ismatrix( model )
    mu = model;
    if size( mu, 1 ) ~= size( y, 1)
        error( 'plotRuleActivation: dim mistmatch mu / y ' )
    end
end

if opts.figure > 0
    h = figure(opts.figure);
else
    h = figure;
end
clf

% time vector given?
if model.ts_ident > 0
    t = [ 0 : model.N_ident-1 ] / model.ts_ident;
else
    t = [ 1 : model.N_ident ];
    opts.labels{1} = 'data-points';
end

% number of diagrams: u and activation (double height)
nu = size(u,2);
ns = nu + 2;
if model.ProductSpace
    % plot y
    ns = ns + 1;
    tit = [opts.title,' in product space z=[u,y]'];
else
    tit = [opts.title,' in input space z=u'];
end

% plot u signals
nf = 0;
for i = 1 : nu
    nf = nf + 1;
    subplot( ns, 1, nf )
    plot( t, u(:,i) )
    grid on
    set(gca, 'XTicklabel', {} )
    if i == 1
        title( tit, 'FontSize',14 )
    end
    if iscell( opts.labels{2} ) == 1
        ylabel( opts.labels{2},'FontSize',14  )
    else
        ylabel( sprintf('%s_{%d}', opts.labels{2},i ),'FontSize',14  )
    end
    
end

if model.ProductSpace
    % plot y signal
    nf = nf + 1;
    subplot(ns,1,nf)
    plot( t, y )
    grid on
    ylabel( opts.labels(3),'FontSize',14  )
    set(gca, 'XTicklabel', {} )
    if ~isempty( opts.limits )
        xlim( opts.limits(1:2) )
    end
end

% plot rule activation
subplot(ns,1,nf+1:nf+2)
hold on
L = {};
for iv= 1 : model.nv
    plot( t, mu(:,iv) )
    L{end+1} = sprintf('v_{%2d}',iv);
end
hold off
grid
box on
set(gca ,'FontSize',14 )
xlabel( opts.labels(1),'FontSize',14  )
ylabel( opts.labels(2),'FontSize',14  )
title( 'Activation \phi_i(z)','FontSize',14  )
legend( L )
% Limits for histogramm (ony x or x&y)
if ~isempty( opts.limits )
    xlim( opts.limits(1:2) )
end
ylim( [-0.05,1.05] )

set(gca ,'FontSize',14 )
xlabel( opts.labels(1),'FontSize',14  )
ylabel( 'activation','FontSize',14  )

h.PaperOrientation = 'landscape';
h.WindowState = 'maximized';
