% TSM_Static_auto

function [ts, yp, par ] = TSM_Static_auto( u,y, par, varargin )

addpath( '../TSModel' )

%% Default limits
% minimal number of clustes / local models
nv_min = 3; 
% minimal data-points per parameter
N_min = 5;  

% Number of data-points
N = size( y, 1 );

%% Default parameters
if nargin < 3 % only u and y -> default par
    
    par.Tries = 10;
    par.nu = size( u, 2 );
    par.nv = nv_min;
    par.ProductSpace = true;
    par.fuzzy = 1.2;
    par.Clustering = 'FCM';
    par.MSF = 'FCM';
    par.LS = 'global';
    par.ParOpt = 'none';  % none/cluster/local/both
    par.IterOpt = 'each'; % each/best
    par.Plots = 'final';
    par.Debug = 1;
    
end

%% Default settings
if ~isfield( par, 'MSF' )
    par.MSF = 'FCM';
end

if ~isfield( par, 'nv' ) || ( isscalar(par.nv) && par.nv == 0 )

    % Minimum x equations per parameter -> nv <= N/x  nPar
    nPar = floor( N / N_min );
    
    if par.ProductSpace
        % nPar = nv*( nu + 1 ) + nv * ( nu + 1 ) = nv * 2 * (nu+1)
        nv_max = max( nv_min, floor( nPar / (2*(par.nu+1)) ) );
    else
        % nPar = nv* nu  + nv * ( nu + 1 ) = nv * 2 * nu
        nv_max = max( nv_min, floor( nPar / (2*par.nu) ) );
    end
    if par.Debug > 0, fprintf('max nv choosen as %d\n', nv_max ),end
    
    par.nv = nv_min : nv_max;
end

% Test mu -> loop mu = [ 1.05, 1.1, 1.2, 1.5 , 2 ]
if ~isfield( par, 'fuzzy' ) || ( isscalar(par.fuzzy) && par.fuzzy == 0 )
    if strcmp( par.MSFt, 'FCM' )
        % FCM
        par.fuzzy = [ 1.05, 1.1, 1.2, 1.5, 2 ];
        if par.Debug > 0, fprintf('choosen range for MSF FCM $\nu$\n' ),end
    elseif strcmp( par.MSFt, 'Gauss' )     
        % Gauss_ sigma
        par.fuzzy = 0;
        if par.Debug > 1, fprintf('choosen range for MSF Gauss $\sigma\n' ),end
    else
        error( 'Unknown MSF type' )
    end
end

if ~isfield( par, 'Clustering' )
    par.Clustering = 'FCM';
end

if ~isfield( par, 'LS' )
    par.LS = 'global';
end

if ~isfield( par, 'Tries' )
    par.Tries = 10;
end

if ~isfield( par, 'ParOpt' )
    par.ParOpt = 'none';
end

if ~isfield( par, 'IterOpt' )
    par.IterOpt = 'best';
end

if ~isfield( par, 'Plots' )
    par.Plots = 'final';
end
if ~isfield( par, 'ProductSpace' )
    par.ProductSpace = true;
end


% Standard labels for inputs and output
for i=1:par.nu
    labels{i} = sprintf('u_{%d}',i );
end
labels{par.nu+1} = 'y';

% Optimize: disable output of lsqnonlin iterations
optimOpt = optimoptions('lsqnonlin');
switch par.Plots
    case 'iter'
        optimOpt.Display = 'iter';
    case 'final'
        optimOpt.Display = 'final';
        optimOpt.Display = 'none';
    otherwise
        optimOpt.Display = 'none';
end

% Best MSE found
best.mse = inf;
best.nv = 0;
best.fuzzy = 0;

%% Loop over nu,m,ber of local models
for nv = par.nv

    
    %% Loop over fuzzy parameter
    for fuzzy = par.fuzzy
       
        if par.Debug > 0
            fprintf( 'Iteration: nv=%2d / fuzzy=%4g\n',nv,fuzzy   )
        end
        tic
        
        %% Inital TS model
        ts = TSModel( 'Static', nv, par.nu, 'comment', 'created by TSM_static' );
        ts.setData( u, y, 'Labels', labels );
        
        nfig = 1;
        for t = 1 : par.Tries
            
            ts.clustering( par.Clustering, 'nue', fuzzy, 'productspace', par.ProductSpace, 'tries', par.Tries );
            
            ts.initialize( par.MSF, 'method', par.LS );
            
            yp = ts.predict( u, y );
            
            
            %% Optimize actual TS model
            if strcmp( par.IterOpt, 'each' )
                
                switch par.ParOpt
                    
                    case 'none'
                    case 'cluster'
                        ts.optimize( 'C', 'optimopts', optimOpt );
                    case 'models'
                        ts.optimize( 'M', 'optimopts', optimOpt );
                    case 'both'
                        ts.optimize( 'B', 'optimopts', optimOpt );
                end
                
                yp = ts.predict( u, y );
                
            end % Opt
            
            dy = y - yp;
            sse = transpose( dy ) * dy;
            mse = sse / length(y);
            if par.Debug > 1
                fprintf( ' try %2d:   mse = %7.4e / best =  %7.4e (nv=%2d/fuzzy=%4g)\n', ...
                    t,mse,best.mse,best.nv,best.fuzzy )
            end
            
            if strcmp( par.Plots, 'iter' )
                % number of model parameters
                np = numel( ts.v ) + numel( ts.A ) + numel( ts.B ) + numel( ts. C );
                ts.plotCluster( ts.v, 'figure', nfig+1,'title', sprintf('Static auto: cluster centers v try %2d',t) );
                plotResiduals( y, yp, 'parameter', np, 'figure', nfig+2, 'title',  sprintf('Static auto: corrleation try %2d',t) );
                nfig = nfig+2;
            end
            
            if mse < best.mse
                best.try = t;
                best.mse = mse;
                best.fuzzy = fuzzy;
                best.nv = nv;
                best.ts = copy( ts );
                best.yp = yp;
            end
            
        end
        
        % Result -> best TS model
        ts = best.ts;
        
        %% Optimize best TS model
        if strcmp( par.IterOpt, 'best' )
            switch par.ParOpt
                
                case 'none'
                case 'cluster'
                    ts.optimize( 'C', 'optimopts', optimOpt );
                case 'models'
                    ts.optimize( 'M', 'optimopts', optimOpt );
                case 'both'
                    ts.optimize( 'B', 'optimopts', optimOpt );
            end
            
            yp = ts.predict( u, y );
            
            dy = y - yp;
            sse = transpose( dy ) * dy;
            mse = sse / length(y);
            best.mse = mse;
            best.yp = yp;
            
            if par.Debug > 1
                fprintf( 'Optimize: mse = %7.4e\n', mse )
            end
            
        end % opt
        if par.Debug > 1, fprintf(' time = %g s\n', toc ), end
         
    end % fuzzy
end % nv

if strcmp( par.Plots, 'final' )
    % number of model parameters
    np = numel( ts.v ) + numel( ts.A ) + numel( ts.B ) + numel( ts. C );
    h.PaperOrientation = 'landscape';
    h = ts.plotCluster( ts.v, 'figure', 1, 'title', 'Static auto: cluster centers v' );
    h.WindowState = 'maximized';
    h = plotResiduals( y, best.yp, 'parameter', np, 'figure', 2, 'title', 'Static auto: correlation' );
    h.PaperOrientation = 'landscape';
    h.WindowState = 'maximized';
end
