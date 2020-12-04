%% PLot MSF in 2D for each dimension as subplots
%%
% Inputs:
%%
%  ts
%  Optional:
%    msf
%    figure
%    split    new figure each split subplots
%%
% Outputs;
%%
%  h    handle to graphics

function h = tsm_plot_msf( ts, varargin )

p = inputParser;
valFcn = @(x) isa(x,'TSModel')
p.addRequired( 'ts', valFcn )
p.addParameter( 'c', [], @ismatrix )
p.addParameter( 'm', 20, @isscalar )
p.addParameter( 'msf', [],@isvector )
p.addParameter( 'limits', [],@isvector )
p.addParameter( 'figure', [], @isscalar )
p.addParameter( 'n', 0, @isscalar )
%p.addParameter( 'split', 0, @isscalar )
p.parse( ts, varargin{:} )
opts = p.Results;

if isempty( opts.c )
    opts.c = ts.c;
else
    if ~isequal( size(opts.c), [ts.nc, ts.nz] )
        error(' tsmplt_msf: c')
    end
    c = opts.c;
end

if isempty( opts.msf )
    % Default: plot all MSF
    nz = ts.nz;
    opts.msf = 1:nz;
else
    nz = size( opts.mfs );
end

% if opts.split > 0
%     np = opts.split;
% end

if opts.n == 0
    opts.n = 10;
end

% Limits u/y -> z_lag_u,z_lag_y
if isempty( opts.limits )
    opts.limits = ts.Limits;
end

if size(opts.limits,1) ~= ts.nc || size(opts.limits,2) < nz 
    error( 'tsm_plt_msf: Limits' )
end

z = makeGrid(  opts.n, nz,opts.limits );
mu = ts.z_msf( z, opts.c, opts.m );

if opts.figure
    h =figure( opts.figure)
else
    h = figure;
end
clf

[sr,sc] = getSubplotPar( nz,1 );
is = 0;
l = {};
for ir=1:sr
    for ic=1:sc
        is = is + 1;
        subplot(sr,sc,is);
        plot( z(:,is), mu(:,is),'-')
        for ic=1:ts.nc
           line(opts.c(ic,is)*[1,1],[0,1],'color','r')
           text(opts.c(ic,is),1,sprintf('c_{%d}',ic) )
        end
        ylim([-0.05,1.05])
        axis square
       
        grid on
        box on
        xlabel( ts.Labels{is} )
        ylabel( 'mu(z)' )
    end
end
