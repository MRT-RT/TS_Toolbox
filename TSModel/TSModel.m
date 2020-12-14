%% Class: TSModel
%
% $Id$ 11.5.2020/ad

% ToDo: Limits u=min/max or Limits lsnonlin
% 29.3.20/[obj.A,obj.B] als Phi?

%%
% Compute and optimize a nonlinear Takagi--Sugeno--Model
% by clustering and identification opf local models
%

%% Given: MISO system
%
% * input vector $u(n x n_u)$ (MI)
% * output vector $y(n x 1)$ (SO)
% * number of clusters $n_c$
% * clustering algoithm (FCM/GK)
% * membership function type (FBF/Gauss) with fuziness parameter $\nue$
% * type of local models (LiP/ARX/OE)
% * lags for scheduling variables $z_{lag}_u}$ and $z_{lag}_y}$ (default = $[1\ldots n_u]$ and $[1\ldots n_y]$
% * lags for regressor variables $x_{lag}_u}$ and $x _{lag}_y}$ (default = $[1\ldots n_u]$

%% Mathematics
%
% $$\hat{y}(k) = \sum_{i=1}^{n_c} \mu_i(z(k))\cdot \hat{y}_{li}(k)$$
%
% with the scheduling variables
%
% $$z = \left[ z_{{lag}_u}, z_{{lag}_y}\right)]$$
%
% and the local models
%
% $$\hat{y}_{li}(k) = \cdot A_i\cdot y(x) + B_i\cdot U(x) + C_i$$
%
% with the regressor
%
% $$x = \left[ x_{lag}_u}), x_{lag}_u} \right]$$

%% Type: LiP-Model
%
% $A=[]$, $z = u(1\ldots n_u)$

%% Type: ARX-Model
%
% $z = u(1\ldots n_u)$

%% Type: OE-Model
%
% Initialize A,B,C randomly or with ARX model


classdef TSModel  < handle & matlab.mixin.Copyable
    
    % ToDo: LiP: kein t -> 1:n
    
    properties
        
        P = struct( ...
            'Version', '1.3', ...
            'Date', '12.11.2020', ...
            'Author', 'Axel Dürrbaum', ...
            'Organisation', 'University of Kassel / ISAC / MRT', ...
            'eMail', 'axeld@uni-kassel.de' )
        
        %% Model order
        Type      % Type of TS model = LiP/NARX/NOE
        nc        % number of clusters
        nu        % number of inputs
        ny = 1    % number of ouputs
        
        %% Data for model identification (Clustering/LS)
        n_ident = 0    % number of data points
        t_ident = [];  % time vector (n x 1 )
        u_ident = []   % vector(n x nu)
        y_ident = []   % vector(n x 1)
        ts_ident = 1;  % sampling time
        Labels = {}    % { u_1,...,u_nu, y }
        Limits = []    % [ u_1,...,u_nu, y x 2 ] lower bound/upper bound
        C_ident = ''   % Comment to dataset
        
        %% Scheduling variable z = (u|y)
        z              % scheduling matrix
        nz             % number of scheduling variables
        z_lag_u = {}   % input lags  0 ... n-1
        z_lag_nu = 0   % number of input lags
        z_lag_y = []   % output lags 1 ... n
        z_maxlag       % maximal lag
        
        z_fct = []     % pointer to sched var function z = fct(u,y)
        z_msf = []     % pointer to membership function msf(u)
        z_Type         % type of  membership function msf(u) Fuzzy/Gauss
        
        zmu            % membership degrees of z
        
        
        %% Cluster c
        c_Type         % type of clustering: fcm/gk/kmeans/...
        c              % vector of initial clusters (nc x nu )
        ProductSpace = false % Clustering in product space [u,y]
        
        Seed = Inf     % Seed fo random number generator
        
        nue            % FBF Fuzziness parameter
        sigma          % Gauss
        m              % FBF: 2 / ( nue - 1 ) / Gauss: -1/(2*sigma^2)
        
        % Parameter for FCM clustering
        FCM_par = struct( 'Exponent', 1, 'MaxIt', 100, ...
            'MinImprove', 1e-5, 'Display', true )
        
        % Parameter for GK clustering
        GK_par = struct( 'Tolerance', '1e-5', 'Display', 'iter' )
        
        % Parameter for KMeans clustering
        KMeans_par = struct( 'Display', 'iter', 'Distance', 'sqeuclidian' )
        
        
        %% Regressor x = [ y(lag_y) | u(lag_u) ]
        x         % regressor matrix
        nx        % number of regressor variables
        x_lag_u = {}  % input lags  0 ... n-1
        x_lag_nu = 0   % number of input lags
        x_lag_y = []  % output lags 1 ... n-1
        x_weights % x <- W * x
        x_maxlag
        
        x_fct = [] % pointer to regressor var function x = fct(u,y)
        
        %% Local model parameter
        A = []     % matrix nc x x_lag_y
        nA         % number of variables len(x_lag_y)
        B = []     % matrix nc x nu * x_lag_u / LS: matrix nc x nu
        nB         % number of variables len(x_lag_u)
        C = []     % matrix nc x 1
        
        Theta = [] % list of local model parameters
        l_Type     % Type of local model initialisation LS local/LS global/random/manual
        
        maxlag = 0 % max( z_maxlag, x_maxlag )
        
        %% Parameter and results of optimization
        OptimPR
        o_Type      % Type of optomizing local model parameters MF/LM/MF+LM
        
        %%
        Name = ''
        Comment = 'TS model'
        Date = now
        
    end
    
    methods
        
        function obj = TSModel( type, nc, nu, varargin )
            
            v = ver('MATLAB');
            if v.Version < 9.7 % R2019b
                error( 'TSModel: need at least Matlab R2019b' )
            end
            
            addpath( '../Functions' )
           
            p = inputParser;
            valFcn = @(x) isscalar(x) && x>0;
            p.addRequired( 'type', @ischar )
            p.addRequired( 'nc', valFcn )
            p.addRequired( 'nu', valFcn )
            p.addParameter('Name','',@ischar)
            p.addParameter('Comment','',@ischar)
            p.addParameter('z_lag_u',[],@ismatrix)
            p.addParameter('z_lag_y',[],@ismatrix)
            p.addParameter('x_lag_u',[],@ismatrix)
            p.addParameter('x_lag_y',[],@ismatrix)
            p.parse( type, nc, nu, varargin{:} )
            opts = p.Results;
            
            obj.Type = type;
            obj.nu = nu;
            obj.nc = nc;
            
            switch type
                
                case 'LiP'
                    obj.z_fct = @tsm_sched_LiP;
                    obj.nz = obj.nu; % MSF only in input-space
                    obj.x_fct = @tsm_reg_LiP;
                    obj.nx = obj.nu + 1; % (B+C)
                    for i=1:obj.nu
                        obj.z_lag_u{i} = [0];
                        obj.x_lag_u{i} = [0];
                    end
                    obj.z_lag_y = [0];
                    obj.x_lag_y = [0];
                    
                case { 'ARX', 'OE' }
                    obj.setSchedulingLags( opts.z_lag_u, opts.z_lag_y);
                    obj.z_fct = @tsm_sched_ARX;
                    obj.setRegressorLags( opts.x_lag_u, opts.x_lag_y);
                    obj.x_fct = @tsm_reg_ARX;
                    obj.nx = length( obj.x_lag_u) +length( obj.x_lag_y); % A+B+C
                
                otherwise
                    error( 'TSModel_ unnknown type <%s>',type )
            end
            obj.Name = opts.Name;
            obj.Comment = opts.Comment;
            obj.Date = datetime('now');
            
        end
        
        %% Set same lags for scheduling and regressor
        function obj = setLags( obj, u_lags, y_lags )
            obj.setSchedulingLags( u_lags, y_lags );
            obj.setRegressorLags( u_lags, y_lags )
        end
                
        function obj = setSchedulingLags( obj, u_lags, y_lags )
            if isempty( u_lags )
                for i=1:obj.nu
                    obj.z_lag_u{i} = 0;
                end
            else
                % Check u is cell array
                % Check dim u = nu x ?
                if size( u_lags,2 ) ~= obj.nu
                    error( 'TSModel/setSchedulingLags: dim lags u' )
                end
                if iscell( u_lags )
                    obj.z_lag_u = u_lags;
                else
                    obj.z_lag_u = {u_lags};
                end
            end
            % Number of input lags in scheduling
            obj.z_lag_nu = 0;
            for i= 1:obj.nu
                obj.z_lag_nu = obj.z_lag_nu + length( obj.z_lag_u{i} );
            end
            if isempty( y_lags )
                obj.z_lag_y = 1:obj.ny;
            else
                % Check dim y = ny x 1?
                if size( y_lags,1 ) ~= 1
                    error( 'TSModel/setSchedulingLags: dim lags y' )
                end
                obj.z_lag_y = y_lags;
            end
            
            obj.z_maxlag = 0;
            obj.maxlag = max( obj.maxlag, obj.z_maxlag );
            obj.nz = length( obj.z_lag_y);
            
            for iu = 1:obj.nu
                obj.z_maxlag = max( obj.z_maxlag, max(obj.z_lag_u{iu}) );
                obj.nz = obj.nz + length( obj.z_lag_u{iu} );
            end
            obj.z_maxlag = max( obj.z_maxlag, max(obj.z_lag_y) );
        end
        

        %% Set regressor lags for u / y
        function obj = setRegressorLags( obj, u_lags, y_lags )
            if isempty( u_lags )
                for i=1:obj.nu
                    obj.x_lag_u{i} = 0;
                end
            else
                % Check u is cell array
                % Check dim u = nu x ?
                if size( u_lags,2 ) ~= obj.nu
                    error( 'TSModel/setRegressorLags: dim lags u' )
                end
                if iscell( u_lags )
                    obj.x_lag_u = u_lags;
                else
                    obj.x_lag_u = {u_lags};
                end
            end
            % Number of input lags in regressor
            obj.x_lag_nu = 0;
            for i= 1:obj.nu
                obj.x_lag_nu = obj.x_lag_nu + length( obj.x_lag_u{i} );
            end
            if isempty( y_lags )
                obj.x_lag_y = 1:obj.ny;
            else
                % Check dim y = ny x ?
                if size( y_lags,1 ) ~= 1
                    error( 'TSModel/setRegressorLags: dim lags y' )
                end
                obj.x_lag_y = y_lags;
            end
            
            obj.nA = size( obj.x_lag_y, 2 );
            obj.x_maxlag = max( obj.x_lag_y);
            obj.maxlag = max( obj.maxlag, obj.x_maxlag );

            obj.nB = 0;
            for iu = 1: obj.nu
                obj.nB = obj.nB + size( obj.x_lag_u{iu}, 2);
                obj.x_maxlag = max( obj.x_maxlag, max(obj.x_lag_u{iu}) );
            end
        end
        
        function obj = setDataLabels( obj, labels )
            if ~isempty( labels )
                if length( labels ) ~= obj.nu + obj.ny
                    error( 'tsModel/Data: Labels')
                end
                obj.Labels = labels;
            end
        end
        
        function obj = setDataComment( obj, comment )
            obj.C_ident = comment;
        end
        
        function obj = setDataLimits( obj, limits )
            if isempty( limits )
                obj.Limits(:,1) = [min(obj.u_ident),min(obj.y_ident)];
                obj.Limits(:,2) = [max(obj.u_ident),max(obj.y_ident)];
            else
                if isequal( size(limits), [obj.nu+obj.ny,2])
                    obj.Limits = limits;
                else
                    error( 'tsModel/Data: Limits')
                end
            end
        end
        
        function obj = setData( obj, u, y, varargin )
            
            p = inputParser;
            p.addRequired( 'u', @ismatrix )
            p.addRequired( 'y', @ismatrix )
            p.addParameter( 't',[],@isvector)
            p.addParameter( 'SampleTime',1,@isscalar)
            p.addParameter( 'Labels',{},@iscell)
            p.addParameter( 'Limits',[],@ismatrix)
            p.addParameter( 'Comment','',@ischar)
            p.parse(u, y, varargin{:} )
            opts = p.Results;
            
            obj.ts_ident = opts.SampleTime;
            
            [n,n_u] = size( u );
            if n_u ~= obj.nu
                error( 'TSModel/setData:  dim u <> n x nu' )
            end
            obj.n_ident = n;
            obj.u_ident = u;
            
            if ~isempty( opts.t )
                [r,c] = size( opts.t );
                if r == 1
                    t = transpose(t);
                end
                if  c == obj.n
                    obj.t_ident = opts.t;
                else
                    error( 'TSModel/Data: t')
                end
            end
            
            [n,n_y] = size( y );
            if n_y ~= 1 || n ~= obj.n_ident
                error( 'TSModel: dim y not n x 1' )
            end
            obj.y_ident = y;
            
            obj = setDataLabels( obj, opts.Labels );
            obj = setDataLimits( obj, opts.Limits );
            obj = setDataComment( obj, opts.Comment );
            
            obj.t_ident = [ 0:obj.n_ident-1 ]*obj.ts_ident;
            
        end
        
        function obj = setCluster( obj, c )
            obj.c = c;
        end
        
        function c = getCluster( obj )
            c = obj.c;
        end
        
        function obj = setFuziness( obj, nue )
            obj.nue = nue;
            obj.m = 2 / (obj.nue-1 );
        end
        
        function obj = clustering( obj, type, varargin )
            
            p = inputParser;
            p.addRequired( 'type', @ischar )
            p.addParameter('nue',0,@isscalar)          % fcm
            p.addParameter('tolerance',1e-5,@isscalar) % gk
            p.addParameter('seed',Inf,@isscalar)
            p.addParameter('tries',1,@isscalar)
            p.addOptional('productspace',false,@islogical)
            
            p.parse( type, varargin{:} )
            opts = p.Results;
            
            % Scheduling variable z_lag_u/z_lag_y or fct()
            obj.c_Type = type;
            obj.ProductSpace = opts.productspace;
            
            obj.nue = opts.nue;
            obj.m = 2 / (obj.nue-1 );
            
            obj.z = obj.z_fct( obj );
            
            if ~isinf( opts.seed )
                obj.Seed = opts.seed;
                rng( opts.seed );
            end
            
            switch opts.type
                case 'FCM'
                    if opts.nue <= 1
                        error( 'tsModel/clustering: nue < 1')
                    end
                    best = inf;
                    fcmopt = [ obj.nue, obj.FCM_par.MaxIt,...
                        obj.FCM_par.MinImprove,obj.FCM_par.Display ];
                    for i=1:opts.tries
                        [ci,~,objFunc ] = fcm( obj.z, obj.nc, fcmopt);
                        if objFunc(end) < best
                            best = objFunc(end);
                            obj.c = ci;
                        end
                    end
                    
                case 'GK'
                    obj.c = gk( obj.z, obj.nc, obj.m, obj.GK_par.Tolerance );
                
                case 'KMeans'
                    [ ~, obj.c ] = kmeans( obj.z, obj.nc, ...
                        'Display', obj.KMeans_par.Display );
                
                otherwise
                    error( 'tsModel/clustering: unknown type <%s>', opts.type)
            end
            
            % Strip y from cluster dimensions
            if obj.ProductSpace
                obj.c = obj.c( :, 1:obj.nu );
            end
        end
        
        %% Initialize model parameters A/B/C
        function obj = initialize( obj, msf, varargin )
            
            p = inputParser;
            p.addRequired( 'msf', @ischar )
            p.addParameter( 'method', 'global', @ischar )
            % FBF
            p.addParameter( 'nue', 1.2, @isscalar)
            % Gauss
            p.addParameter( 'sigma', 2, @isscalar)
            
            p.parse( msf, varargin{:} )
            opts = p.Results;
            
            switch msf
                case 'FBF'
                    obj.z_Type = 'F';
                    obj.nue = opts.nue;
                    obj.m = 2 / ( obj.nue-1 );
                    obj.z_msf = @tsm_membership_FBF;
                case 'Gauss'
                    obj.z_Type = 'G';
                    obj.sigma = opts.sigma;
                    obj.m = -1/(2*opts.sigma^2);
                    obj.z_msf = @tsm_membership_Gauss;
                otherwise
                    error( 'TSModel/initialize: unknown MSF <%s>', opts.msf)
            end
            
            switch obj.Type
                case 'LiP'
                    
                    % Scheduling var = u
                    obj.z = obj.z_fct( obj );
                    obj.zmu = obj.z_msf( obj.z, obj.c, obj.m );
                    obj.z = [ obj.z, ones(obj.n_ident,1) ];
                    
                    obj.l_Type = opts.method;
                    %Todo: random/manual
                    
                    if strcmp( opts.method, 'local' ) % local LS
                        obj.Theta = [];
                        for ic = 1 : obj.nc
                            Phi = bsxfun( @times, obj.z, obj.zmu(:,ic)  );
                            theta = Phi \ obj.y_ident;
                            obj.Theta = [obj.Theta, theta ];
                            phi =  transpose( reshape(  theta, obj.nu+1,1 ) );
                            obj.B( ic, : ) = phi( 1:obj.nu );
                            obj.C( ic, 1 ) = phi( obj.nu+1 );
                        end
                        
                    else % global LS
                        
                        Phi = bsxfun(@times, obj.z, obj.zmu(:,1) );
                        for ic = 2 : obj.nc
                            Phi = [ Phi, bsxfun(@times, obj.z, obj.zmu(:,ic) ) ];
                        end
                        obj.Theta = Phi \ obj.y_ident;
                        phi = transpose(  reshape( obj.Theta,obj.nu+1,obj.nc) );
                        obj.B = phi( :, 1:obj.nu );
                        obj.C = phi( :, obj.nu+1 );
                        
                    end
                    
                case { 'ARX', 'OE' }
                    
                    % Scheduling var z = [ y(lags_y) | u(lags_u) | 1 ]
                    obj.z = obj.z_fct( obj );
                    n = size( obj.z,1);
                    obj.zmu = obj.z_msf( obj.z, obj.c, obj.m );
                    obj.z = [ obj.z, ones(n,1) ];

                    %ToDo: #z_lags <> #x_lags???
                    if  strcmp( opts.method, 'global' ) % local LS
                        
                        Phi = bsxfun(@times, obj.z, obj.zmu(:,1) );
                        for ic = 2 : obj.nc
                            Phi = [ Phi, bsxfun(@times, obj.z, obj.zmu(:,ic) ) ];
                        end
                        obj.Theta = Phi \ obj.y_ident(obj.z_maxlag+1:end);
                        phi = transpose(  reshape( obj.Theta ,obj.nA+obj.nB+1,obj.nc) );
                        obj.A = phi( :, 1:obj.nA );
                        obj.B = phi( :, obj.nA+(1:obj.nB) );
                        obj.C = phi( :, obj.nA+obj.nB+1 );

                    elseif  strcmp( opts.method, 'local' )

                        for ic = 1 : obj.nc
                            obj.Phi = bsxfun( @times, obj.z, obj.zmu(:,ic)  );
                            phi =  transpose( reshape(  obj.Phi \ obj.y_ident(1:n), obj.nA+obj.nB+1,1 ) );
                            obj.A(ic,:) = phi( 1:obj.nA );
                            obj.B(ic,:) = phi( obj.nA+(1:obj.nB) );
                            obj.C(ic,:) = phi( obj.nA+obj.nB+1 );
                        end
                        
                    else
                        error('TSModel/initialize: unknowm LS method <%s>', obj.method)
                    end
            end
        end
        
        function [A,B,C] = getLM( obj )
            A = obj.A;
            B = obj.B;
            C = obj.C;
        end
        
        function obj = setLM( obj, A,B,C )
            
            if ~isempty(A) && ~isequal( size(A), [obj.nc,length(obj.x_lag_y)] )
                error( 'TSModel/setLM: dim A <> nc x lag_y' )
            end
            obj.A = A;
            obj.nA = size(A,2);
            
            if ~isequal( size(B), [obj.nc,length(obj.x_lag_u)] )
                error( 'TSModel/setLM: dim B <> nu x nc x lag_u' )
            end
            
            obj.B = B;
            obj.nB = size(B,2);
            if ~isequal( size(C), [obj.nc,1] )
                error( 'TSModel/setLM: dim C <> nc x 1' )
            end
            obj.C = C;
        end
        
        function obj = set_msf_type( obj, type )
            switch type
                case 'FBF'
                    obj.z_msf = @tsm_membership_FBF;
                case 'Gauss'
                    obj.z_msf = @tsm_membership_Gauss;
                otherwise
                    error( 'tsModel/set_msf: unknown type <%s>',type)
            end
        end
        
        %% Evaluate system at vectors u / u,y
        function yp = evaluate( obj, u, y )

            if isempty(obj.B)
                error( 'tsModel/evaluate: matrix B empty' )
            end
            if isempty(obj.C)
                error( 'tsModel/evaluate: matrix C empty' )
            end
            switch obj.Type
                
                case 'LiP'
                    % Check A/B/C not empty (model not yet initialized)
                    yp = tsm_evaluate_LiP( obj, u );
                case 'ARX'
                    if isempty(obj.A) || isempty(obj.C)
                        error( 'tsModel/evaluate: matrix A empty' )
                    end
                    yp = tsm_evaluate_ARX( obj, u, y );
                case 'OE'
                    yp = tsm_evaluate_OE( obj, u, y );
            end
        end
        
        %%
        function obj = optimize( obj, method, varargin )
            
            p = inputParser;
            p.addRequired( 'method', @ischar )
            p.addParameter('optimopts',[],@isstruct)
            p.parse( method, varargin{:} )
            opts = p.Results;
            
            % Options for LSQNONLIN
            if isempty( opts.optimopts )
                optimopts = optimoptions('lsqnonlin');
                optimopts.Display = 'Iter';
                %optimopts.Algorithm = 'levenberg-marquardt';
                optimopts.Algorithm = 'trust-region-reflective';
                optimopts.FunctionTolerance = 1e-3;
                optim.OptimalityTolerance = 1e-3;
                optimopts.StepTolerance = 1e-6;
            else
                optimopts = opts.optimopts;
            end
            
            switch obj.Type
                
                case 'LiP'
                    switch method
                        case 'M' % optimize membership MF
                            % Limit c to range of u|y
                            obj.o_Type = 'MF';
                            lb = reshape( repmat(obj.Limits(1:obj.nu,1),1,obj.nc),1,obj.nc*obj.nu );
                            ub = reshape( repmat(obj.Limits(1:obj.nu,2),1,obj.nc),1,obj.nc*obj.nu );
                            x0 = [ reshape(obj.c,1,obj.nc*obj.nu) ];
                            [x_opt,resnorm,residual,exitflag,output] = ...
                                lsqnonlin( @(x) tsm_optimize_LiP_MF( x,...
                                obj.u_ident, obj )-obj.y_ident,...
                                x0, lb,ub,optimopts);
                            obj.c = reshape( x_opt(1:obj.nc*obj.nu), obj.nc, obj.nu );
                        case 'L' % optimize local models LM
                            obj.o_Type = 'LM';
                            lb = [];
                            ub = [];
                            x0 = [ reshape(obj.B,1,obj.nc*obj.nu),...
                                reshape(obj.C,1,obj.nc) ];
                            [x_opt,resnorm,residual,exitflag,output] = ...
                                lsqnonlin( @(x) tsm_optimize_LiP_LM( x,...
                                obj.u_ident, obj )-obj.y_ident,...
                                x0, lb,ub,optimopts);
                            n1 = 1; n2 = obj.nc * obj.nu;
                            obj.B = reshape( x_opt(n1:n2), obj.nc,obj.nu );
                            n1 = n2 + 1; n2 = n1+obj.nc-1;
                            obj.C = reshape( x_opt(n1:n2),obj.nc,1 );
                        case 'B' % optimize MF and LM
                            obj.o_Type = 'MF&LM';
                            lb = [ reshape( repmat(obj.Limits(1:obj.nu,1),1,obj.nc),1,obj.nc*obj.nu ),...
                                -inf*ones(1,obj.nc*(obj.nu+1)) ];
                            ub = [ reshape( repmat(obj.Limits(1:obj.nu,2),1,obj.nc),1,obj.nc*obj.nu ),...
                                +inf*ones(1,obj.nc*(obj.nu+1)) ];
                            x0 = [ reshape(obj.c,1,obj.nc*obj.nu),...
                                reshape(obj.B,1,obj.nc*obj.nu),...
                                reshape(obj.C,1,obj.nc) ];
                            [x_opt,resnorm,residual,exitflag,output] = ...
                                lsqnonlin( @(x) tsm_optimize_LiP( x,...
                                obj.u_ident, obj )-obj.y_ident,...
                                x0, lb,ub,optimopts);
                            n2 = obj.nc * obj.nu;
                            obj.c = reshape( x_opt(1:n2), obj.nc, obj.nu );
                            n1 = n2 + 1; n2 = 2*obj.nc*obj.nu;
                            obj.B = reshape( x_opt(n1:n2), obj.nc,obj.nu );
                            n1 = n2 + 1; n2 = n1+obj.nc-1;
                            obj.C = reshape( x_opt(n1:n2),obj.nc,1 );
                        otherwise
                            error('TSModel/optimize: not M, L, or B' )
                    end
                    obj.OptimPR.optimopts = optimopts;
                    obj.OptimPR.exitflag = exitflag;
                    obj.OptimPR.resnorm = resnorm;
                    obj.OptimPR.resitual = residual;
                    obj.OptimPR.output = output;
                    
                case 'ARX'
                    switch method
                        case 'M' % optimize membership MF
                            % Limit c to range of u|y
                            obj.o_Type = 'MF';
                            nu = obj.nB * obj.nu; % Lags beachten
                            lb = reshape( repmat(obj.Limits(1:obj.nu,1),1,obj.nc),...
                                1,obj.nc*obj.nu );
                            ub = reshape( repmat(obj.Limits(1:nu,2),1,obj.nc),1,obj.nc*obj.nu );
                            x0 = [ reshape(obj.c,1,obj.nc*obj.nu) ];
                            [x_opt,resnorm,residual,exitflag,output] = ...
                                lsqnonlin( @(x) tsm_optimize_ARX_MF( x,...
                                obj.u_ident, obj )-obj.y_ident,...
                                x0, lb,ub,optimopts);
                            obj.c = reshape( x_opt(1:obj.nc*obj.nu), obj.nc, obj.nu );
                        case 'L' % optimize local models LM
                            obj.o_Type = 'LM';
                            lb = [];
                            ub = [];
                            x0 = [ reshape(obj.B,1,obj.nc*obj.nu),...
                                reshape(obj.C,1,obj.nc) ];
                            [x_opt,resnorm,residual,exitflag,output] = ...
                                lsqnonlin( @(x) tsm_optimize_ARX_LM( x,...
                                obj.u_ident, obj )-obj.y_ident,...
                                x0, lb,ub,optimopts);
                            n1 = 1; n2 = obj.nc * obj.nu;
                            obj.B = reshape( x_opt(n1:n2), obj.nc,obj.nu );
                            n1 = n2 + 1; n2 = n1+obj.nc-1;
                            obj.C = reshape( x_opt(n1:n2),obj.nc,1 );
                        case 'B' % optimize MF and LM ( c | A | B C ]
                            obj.o_Type = 'MF&LM';
                            % lag_y | lag_u
                            lb = [ reshape( repmat(obj.Limits(1:obj.nu,1),1,obj.nc),1,obj.nc*obj.nu ),...
                                -inf*ones(1,obj.nc*(obj.nA+obj.nB+1)) ];
                            ub = [ reshape( repmat(obj.Limits(1:obj.nu,2),1,obj.nc),1,obj.nc*obj.nu ),...
                                +inf*ones(1,obj.nc*(obj.nA+obj.nB+1)) ];
                            
                            p0 = [ reshape(obj.c,1,obj.nc*(obj.nA+obj.nB)),...
                                   transpose(obj.Theta) ];
                           
%                             [x_opt,resnorm,residual,exitflag,output] = ...
%                                 lsqnonlin( @(x) tsm_optimize_ARX( x,...
%                                 obj.u_ident, obj )-obj.y_ident,...
%                                 x0, lb,ub,optimopts);
                            z = tsm_Regressor( obj.u_ident, obj.z_lag_u, obj.y_ident, obj.z_lag_y );
                            y = obj.y_ident(obj.x_maxlag+1:end);
                            [p_opt,resnorm,residual,exitflag,output] = ...
                                lsqnonlin( @(p) tsm_optimize_ARX( p, z, obj )- y,...
                                p0,lb,ub,optimopts );
                            
                            n2 = obj.nc * (obj.nA+obj.nB);
                            obj.c = reshape( p_opt(1:n2), obj.nc, (obj.nA+obj.nB) );
                            obj.Theta = transpose( p_opt(n2+1:end) );

                            n1 = n2 + 1; n2 = n1 + obj.nc*obj.nA-1;
                            obj.A = reshape( p_opt(n1:n2), obj.nc,obj.nA );
                            n1 = n2 + 1; n2 = n1 + obj.nc*obj.nB-1;
                            obj.B = reshape( p_opt(n1:n2), obj.nc,obj.nB );
                            n1 = n2 + 1; n2 = n1+obj.nc-1;
                            obj.C = reshape( p_opt(n1:n2),obj.nc,1 );
                        
                        otherwise
                            error('TSModel/optimize: not M, L, or B' )
                    end
                    obj.OptimPR.optimopts = optimopts;
                    obj.OptimPR.exitflag = exitflag;
                    obj.OptimPR.resnorm = resnorm;
                    obj.OptimPR.residual = residual;
                    obj.OptimPR.output = output;
                    
                case 'OE'
                       obj.o_Type = 'MF&LM';
                            % c = [lag_y,lag_u)
                            %ToDo: u x nu
                             lb = [ repmat( ...
                                    [ repmat(obj.Limits(2,1),1,obj.nA),...   % lag_y
                                      repmat(obj.Limits(1,1),1,obj.nB) ],...  % nu * lag_u
                                      1,obj.nc),... 
                                   -inf*ones(1,obj.nc*(obj.nA+obj.nB+1)) ];
                             ub = [ repmat( ...
                                    [ repmat(obj.Limits(2,2),1,obj.nA),...    % lag_y
                                      repmat(obj.Limits(1,2),1,obj.nB) ],...  % nu * lag_u
                                      1,obj.nc),... 
                                   +inf*ones(1,obj.nc*(obj.nA+obj.nB+1)) ];
                            
                            p0 = [ reshape(obj.c,1,obj.nc*(obj.nA+obj.nB)),...
                                   transpose(obj.Theta) ];
                            
                            [p_opt,resnorm,residual,exitflag,output] = ...
                                lsqnonlin( @(p) tsm_optimize_OE( p, obj ),...
                                p0, lb,ub,optimopts );
                            n2 = obj.nc * (obj.nA+obj.nB);
                            obj.c = reshape( p_opt(1:n2), obj.nc, (obj.nA+obj.nB) );
                            obj.Theta = transpose( p_opt(n2+1:end) );

                            n1 = n2 + 1; n2 = n1 + obj.nc*obj.nA-1;
                            obj.A = reshape( p_opt(n1:n2), obj.nc,obj.nA );
                            n1 = n2 + 1; n2 = n1 + obj.nc*obj.nB-1;
                            obj.B = reshape( p_opt(n1:n2), obj.nc,obj.nB );
                            n1 = n2 + 1; n2 = n1+obj.nc-1;
                            obj.C = reshape( p_opt(n1:n2),obj.nc,1 );
                        
            end
        end
        
        function obj = setName( obj, name )
            obj.Name = name;
        end
        
        function obj = setComment( obj, comment )
            obj.Comment = comment;
        end
        
        function [A,B,C] = Regressor2ABC( obj, Theta )
            if nargin < 2
                Theta = obj.Theta;
            end
            n2 = obj.nc * obj.nA;
            A = reshape(Theta(1:n2), obj.nc, obj.nA);
            n1 = n2+1; n2 = obj.nc * obj.nB;
            B = reshape(Theta(n1:n2), obj.nc, obj.nB);
            n1 = n2+1; n2 = obj.nc * obj.nC;
            C = reshape(Theta(n1:n2), obj.nc, 1);
        end
        
        function disp( obj )
            fprintf( 'TS-Model: Type=%s', obj.Type )
            if ~isempty( obj.Name )
                fprintf( ', Name="%s"', obj.Name )
            end
            if ~isempty( obj.Comment )
                fprintf( ' (%s)', obj.Comment )
            end
            if ~isempty( obj.Date )
                fprintf( ', [created: %s]', obj.Date )
            end
            fprintf( '\n Order: nc = %d, nu = %d, ny = %d\n', obj.nc, obj.nu, obj.ny )
            
            if obj.n_ident > 0
                fprintf( '\n Ident data: n=%d, ts=%g', obj.n_ident, obj.ts_ident )
                if obj.C_ident
                    fprintf( ' (%s)', obj.C_ident )
                end
                fprintf( '\n' )
            end
            fprintf( 'Scheduling lags: ')
            if ~isempty( obj.z_lag_u )
                for i=1:obj.nu
                    fprintf( 'u_%d:%s, ', i,mat2str( obj.z_lag_u{i} ) )
                end
            end
            fprintf( 'y = %s\n', mat2str(obj.z_lag_y) )
            
            fprintf( 'Regressor lags:  ')
            if ~isempty( obj.x_lag_u )
                for i=1:obj.nu
                    fprintf( 'u_%d:%s, ', i,mat2str( obj.x_lag_u{i} ) )
                end
            end
            fprintf( 'y = %s\n', mat2str(obj.x_lag_y) )
            if obj.c_Type
                fprintf( ' Clustering: %s, nue=%g\n', obj.c_Type, obj.nue )
            end
            if obj.l_Type
                fprintf( ' Initialisation local models: LiP %s\n', obj.l_Type )
            end
            if obj.o_Type
                fprintf( ' Optimization parameter: %s\n', obj.o_Type )
            end
        end
        
        %%
        function h = plotCluster( obj, c, varargin )
            
            p = inputParser;
            p.addRequired('c',@ismatrix)
            p.addParameter('figure',2,@isscalar)
            p.addParameter('title','Cluster',@ischar)
            p.addParameter('file','',@ischar)
            p.parse( c, varargin{:} )
            opts = p.Results;
            
            if nargin < 2
                c = obj.c;
            end
            n = size( c, 2 );
            
            h = figure(opts.figure);
            clf
            
            [sr,sc] = getSubplotPar( n );
            is = 0;
            l = {};
            for i1=1:n-1
                for i2=i1+1:n
                    is = is+1;
                    subplot(sr,sc,is)
                    if ~isempty( obj.z )
                        plot( obj.z(:,i1),obj.z(:,i2),'k.' )
                        hold on
                        l{end+1} = 'data';
                    end
                    plot( c(:,i1), c(:,i2), 'rx' )
                    if ~isempty( obj.Limits )
                        axis( [obj.Limits(i1,:), obj.Limits(i2,:) ] )
                    end
                    axis square
                    grid on
                    box on
                    xlabel( obj.Labels(i1) )
                    ylabel( obj.Labels(i2) )
                    if is == 1 && ~isempty( opts.title )
                        title( opts.title )
                        l{end+1} = 'clusters';
                        legend( l, 'location','SE' )
                    end
                end
            end
            
            if ~isempty( opts.file )
                i = strfind( opts.file, '.');
                if i > 0
                    type = opts.file(i+1:end);
                else
                    type = 'png';
                end
                print( '-r200', ['-d',type], opts.file )
            end
        end
        
        function h = plotIdentData( obj, t,u, y  )
            
            %ToDo: u separat/zusammen
            %varargin
            h = figure(1);
            clf
            
            if nargin < 4
                t = obj.t_ident;
                u = obj.u_ident;
                y = obj.y_ident;
            end
            
            n = size(t,2);
            nu = size(u,2);
            ny = size(y,2);
            
            np = obj.nu+obj.ny;
            is = 0;
            for i=1:nu
                is = is+1;
                subplot(nu+ny,1,is)
                plot( t, u(:,i),'k.:' )
                if ~isempty( obj.Limits )
                    ylim( obj.Limits( i,: ) )
                end
                grid on
                ylabel( sprintf('u_{%d}',i ) )
                if is ==1
                    title( sprintf( 'Identification data n=%d', obj.n_ident ) )
                end
            end
            for i=1:ny
                is = is+1;
                subplot(nu+ny,1,is)
                plot( t, y(:,i),'k.:' )
                if ~isempty( obj.Limits )
                    ylim( obj.Limits( i+obj.nu,: ) )
                end
                grid on
                ylabel( sprintf('y_{%d}',i ) )
            end
            xlabel( 'time' )
        end
        
        function plot( obj )
            
            plotIdentData( obj, obj.t_ident, obj.u_ident,  obj.y_ident );
            plotCluster( obj, obj.c );
            
        end
        
    end
end