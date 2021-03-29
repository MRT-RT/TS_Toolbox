%% Class: TSModel
%
% $Id$ 

%%
% Compute and optimize a nonlinear Takagi--Sugeno--Model
% by clustering and identification of local models
%

%% Given: MISO system
%%
%
% * input vector $u(n x n_u)$ (MI)
% * output vector $y(n x 1)$ (SO)
% * number of clusters / local models $n_c$
% * clustering algorithm (FCM/GK)
% * clustering norm (Euclidean/Mahalanobis)
% * membership function type (FCMF/Gauss) with fuzziness parameter $\nu$
% * type of local models (Static/ARX/OE)
% * lags for scheduling variables $z_{{lag}_u}$ and $z_{{lag}_y}$ (default = $[1\ldots n_u]$ and $[1\ldots n_y]$
% * lags for regressor variables $x_{{lag}_u}$ and $x_{{lag}_y}$ (default = $[1\ldots n_u]$

%% Mathematics
% 
% $$\hat{y}(k) = \sum_{i=1}^{n_c} \mu_i(z(k))\cdot \hat{y}_{i}(k)$$
%
% with the scheduling variables
%
% $$z = \left[ y(z_{{lag}_y}), u(z_{{lag}_u}) \right]$$
%
% and the local models
%
% $$ y_{i}(k) = \cdot A_i\cdot y(x_{{lag}_y}) + B_i\cdot u(x_{{lag}_u}) + C_i $$
%
% with the regressor
%
% $$ x = \left[ y(x_{{lag}_y})), u(x_{{lag}_u}) \right] $$

%% Local model type: Static
%
% $A=[]$, $z = u(1\ldots n_u)$

%% Local model type: ARX
%
% $z = u(1\ldots n_u)$

%% Local model type: OE
%
% Initialize A,B,C randomly or with ARX model


classdef TSModel  < tsm_Base
    % TSModel TS model class
    
    properties
        
        P = struct( ...
            'Version', '$Revision$', ...
            'Date', '$Date$', ...
            'Author', 'Axel Dürrbaum', ...
            'Organisation', 'University of Kassel / ISAC / MRT', ...
            'eMail', 'axeld@uni-kassel.de' )
        
        %% Model order
        nv        % number of clusters v
        nu        % number of inputs
        ny = 1    % number of ouputs
        
        %% Premise and conclusion parts
        prem      % premise (tsm_Prem object, clustering type)
        conc      % conclusion (tsm_Conc object, local models type)
        
        %% Data for model identification (Clustering/LS)
        N_ident = 0    % number of data-points
        t_ident = [];  % time vector (n x 1 )
        u_ident = []   % vector(n x nu)
        y_ident = []   % vector(n x 1)
        ts_ident = -1;  % sampling time
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
        z_Type ='FCM'  % type of membership function msf(u) FCM/Gauss
        
        zmu            % membership degrees of z
        
        
        %% Clustering v
        v_Type = 'FCM'       % type of clustering: fcm/gk/kmeans/...
        v_Norm = 'Euclidean' % Norm for clustering (Euclidean/Mahlanobis')
        v              % vector of initial clusters (nv x nu )
        ProductSpace = false % Clustering in product space [u,y]

        
        Seed = Inf     % Seed fo random number generator
        
        nue            % FCM Fuzziness parameter
        sigma          % Gauss
        m              % FCM: 2 / ( nue - 1 ) / Gauss: -1/(2*sigma^2)
        
        % Parameter for FCM clustering
        FCM_par = struct( 'Exponent', 1, 'MaxIt', 100, ...
            'MinImprove', 1e-5, 'Display', false )
        
        % Parameter for GK clustering
        GK_par = struct( 'Tolerance', '1e-5', 'Display', 'iter' )
        
        % Parameter for KMeans clustering
        KMeans_par = struct( 'Display', 'iter', 'Distance', 'qeuclidian' )
        
        
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
        A = []     % matrix nv x x_lag_y
        nA         % number of variables len(x_lag_y)
        B = []     % matrix nv x nu * x_lag_u / LS: matrix nv x nu
        nB         % number of variables len(x_lag_u)
        C = []     % matrix nv x 1
        
        Theta = [] % list of local model parameters
        LS_type    % Type of local model initialisation LS local/LS global/random/manual
        
        maxlag = 0 % max( z_maxlag, x_maxlag )
        
        %% Parameter and results of optimization
        OptimPR
        Opt_Type      % Type of optimizing local model parameters MF/LM/MF+LM
        
    end
    
    methods
        
        function obj = TSModel( type, nv, nu, varargin )
            % constructor 
            
            v = ver('MATLAB');
            if v.Version < 9.7 % R2019b
                error( 'TSModel: need at least Matlab R2019b' )
            end
            
            addpath( '../Functions' )
            
            p = inputParser;
            valFcn = @(x) isscalar(x) && x>0;
            p.addRequired( 'type', @ischar )
            p.addRequired( 'nv', valFcn )
            p.addRequired( 'nu', valFcn )
            p.addParameter('Name','undefined',@ischar)
            p.addParameter('Comment','',@ischar)
            p.addParameter('z_lag_u',{},@iscell)
            p.addParameter('z_lag_y',[],@ismatrix)
            p.addParameter('x_lag_u',{},@iscell)
            p.addParameter('x_lag_y',[],@ismatrix)
            p.parse( type, nv, nu, varargin{:} )
            opts = p.Results;

            obj@tsm_Base( opts.Name, opts.Comment );

            obj.Type = type;
            obj.nu = nu;
            obj.nv = nv;
            
            switch type
                
                case 'Static'
                    obj.z_fct = @tsm_sched_Static;
                    obj.nz = obj.nu; % MSF only in input-space
                    obj.x_fct = @tsm_reg_Static;
                    obj.nx = obj.nu + 1; % (B+C)
                    obj.z_lag_u = {};
                    obj.z_lag_u = [];
                    obj.z_lag_y = {};
                    obj.x_lag_y = [];
                    
                case { 'ARX', 'OE' }
                    obj.setSchedulingLags( opts.z_lag_u, opts.z_lag_y);
                    obj.z_fct = @tsm_sched_ARX;
                    obj.setRegressorLags( opts.x_lag_u, opts.x_lag_y);
                    obj.x_fct = @tsm_reg_ARX;
                    obj.nx = length( obj.x_lag_u) +length( obj.x_lag_y); % A+B+C
                
                otherwise
                    error( 'TSModel_ unnknown type <%s>',type )
            end
                        
        end
        
        function obj = setLags( obj, u_lags, y_lags )
            % Set same lags for scheduling and regressor
            obj.setSchedulingLags( u_lags, y_lags );
            obj.setRegressorLags( u_lags, y_lags );
        end
                
        function obj = setSchedulingLags( obj, u_lags, y_lags )
            % Set lags for scheduling 
            
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

        function obj = setRegressorLags( obj, u_lags, y_lags )
            % Set lags for regressor
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
            % Set comment for identification data
            obj.C_ident = comment;
        end
        
        function obj = setDataLimits( obj, limits )
            % Set limits for identification data
            % limits: [u_1min, u1_max ; ... ; u_nu_min, u_nu_max ; y_min, y_max ]
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
            % Set identification data u and y
            
            p = inputParser;
            p.addRequired( 'u', @ismatrix )
            p.addRequired( 'y', @ismatrix )
            p.addParameter( 't',[],@isvector)
            p.addParameter( 'SampleTime',-1,@isscalar)
            p.addParameter( 'Labels',{},@iscell)
            p.addParameter( 'Limits',[],@ismatrix)
            p.addParameter( 'Comment','',@ischar)
            p.parse(u, y, varargin{:} )
            opts = p.Results;
            
            obj.ts_ident = opts.SampleTime;
            
            [n,n_u] = size( u );
            if n_u ~= obj.nu
                error( 'TSModel/setData: dim u not N x nu' )
            end
            obj.N_ident = n;
            obj.u_ident = u;
            
            if ~isempty( opts.t )
                [r,c] = size( opts.t );
                if r == 1
                    opts.t = transpose( opts.t );
                end
                if  c == obj.N_ident
                    obj.t_ident = opts.t;
                else
                    error( 'TSModel/setData: dim t not N x 1')
                end
            end
            
            [n,n_y] = size( y );
            if n_y ~= 1 || n ~= obj.N_ident
                error( 'TSModel/setData: dim y not N x 1' )
            end
            obj.y_ident = y;
            
            obj = setDataLabels( obj, opts.Labels );
            obj = setDataLimits( obj, opts.Limits );
            obj = setDataComment( obj, opts.Comment );
            
            % ts_ident vs. t_ident
            if opts.SampleTime > 0 
                % create missing time vector from sample-time
                if isempty( obj.t_ident )
                    obj.t_ident = [ 0:obj.N_ident-1 ]*obj.ts_ident;
                else
                    if opts.SampleTime ~= diff(obj.t_ident(1:2))
                        error( 'TSModel/setData: sample-time mismatch' )
                    end
                end
            end
            
        end
        
        function obj = setCluster( obj, arg )
            % setCluster
            % Input: arg = nv (Scalar) or v (matrix)
            if isscalar( arg )
                obj.nv = arg;
                obj.v = zeros( obj.nv, obj.nu )
            elseif ismatrix( arg )
                obj.v = arg;
                [obj.nv,nu] = size(arg,2);
                % #columns / ProductSpace
                if (obj.ProductSpace && nu ~= obj.nu + 1) || (nu ~= obj.nu)
                    error( 'tsModel/setCluster: dim error v col<>nu' )
                end
            else
                error( 'tsModel/setCluster: arg not nv or v' )
            end
        end
        
        function v = getCluster( obj )
            % Get cluster centers v
            v = obj.v;
        end
        
        function obj = setFuzziness( obj, nue )
            % Set fuzziness factor (FBF or Gauss)
            obj.nue = nue;
            obj.m = 2 / (obj.nue-1 );
        end
        
        function obj = clustering( obj, type, varargin )
            % Do the clustering
            p = inputParser;
            p.addRequired( 'type', @ischar )
            p.addParameter('norm',obj.v_Norm,@ischar)  % fcm norm
            p.addParameter('tries',1,@isscalar)
            p.addParameter('nue',obj.nue,@isscalar)    % fcm
            p.addParameter('tolerance',1e-5,@isscalar) % gk
            p.addParameter('seed',Inf,@isscalar)
            p.addOptional('productspace',false,@islogical)
            
            p.parse( type, varargin{:} )
            opts = p.Results;
            
            % Scheduling variable z_lag_u/z_lag_y or fct()
            obj.v_Type = opts.type;
            obj.v_Norm = opts.norm;
            obj.ProductSpace = opts.productspace;
            
            obj.nue = opts.nue;
            obj.m = 2 / (obj.nue-1 );
            
            obj.z = obj.z_fct( obj );
            
            if ~isinf( opts.seed )
                obj.Seed = opts.seed;
                rng( opts.seed );
            end
            
            switch obj.v_Type
                case 'FCM'
                    if opts.nue <= 1
                        error( 'tsModel/clustering: nue < 1')
                    end
                    best = inf;
                    fcmopt = [ obj.nue, obj.FCM_par.MaxIt,...
                        obj.FCM_par.MinImprove,obj.FCM_par.Display ];
                    switch obj.v_Norm
                        case 'Euclidean'
                            for i=1:opts.tries
                                [vi,~,objFunc ] = fcm_Euclidean( obj.z, obj.nv, fcmopt);
                                if objFunc(end) < best
                                    best = objFunc(end);
                                    obj.v = vi;
                                end
                            end
                        case 'Mahalanobis'
                            for i=1:opts.tries
                                [vi,~,objFunc ] = fcm_Mahalanobis( obj.z, obj.nv, fcmopt);
                                if objFunc(end) < best
                                    best = objFunc(end);
                                    obj.v = vi;
                                end
                            end
                    end
                    
                case 'GK' % Gustaffson-Kessel
                    obj.v = gk( obj.z, obj.nv, obj.m, obj.GK_par.Tolerance );
                
                case 'KMeans'
                    [ ~, obj.v ] = kmeans( obj.z, obj.nv, ...
                        'Display', obj.KMeans_par.Display );
                
                otherwise
                    error( 'tsModel/clustering: unknown type <%s>', opts.type)
            end
            
            % Strip y from cluster/sched dimensions
            if obj.ProductSpace
                obj.v = obj.v( :, 1:obj.nu );
                obj.z = obj.z( :, 1:obj.nu );
            end
        end
        
        function mu = getMSF( obj, u, y )
            % Get membership function degree
            
            if nargin < 2
                u = obj.u_ident;
                y = obj.y_ident;
            end
            
            switch obj.Type
                case 'Static'
                    % Scheduling var = u / forget ProductSpace
                    obj.z = obj.z_fct( obj, u, y );
                    mu = obj.z_msf( obj.z(:,1:obj.nu), obj.v, obj.m ); % Skip y ProdoctSpace
                case { 'ARX', 'OE' }
                    % Scheduling var z = [ y(lags_y) | u(lags_u) | 1 ]
                    obj.z = obj.z_fct( obj, u );
                    n = size( obj.z,1);
                    mu = obj.z_msf( obj.z, obj.v, obj.m );
            end
        end
        
        function obj = initialize( obj, msf, varargin )
            % Initialize local models: parameters A/B/C 
            % with local or global Least Squares
            
            p = inputParser;
            p.addRequired( 'msf', @ischar )
            p.addParameter( 'method', 'global', @ischar )
            % FCM
            p.addParameter( 'nue', 1.2, @isscalar)
            % Gauss
            p.addParameter( 'sigma', 2, @isscalar)
            
            p.parse( msf, varargin{:} )
            opts = p.Results;
            
            switch msf
                case 'FCM'
                    obj.z_Type = 'FCM';
                    obj.nue = opts.nue;
                    obj.m = 2 / ( obj.nue-1 );
                    obj.z_msf = @tsm_membership_FCM;
                case 'Gauss'
                    obj.z_Type = 'Gauss';
                    obj.sigma = opts.sigma;
                    obj.m = -1/(2*opts.sigma^2);
                    obj.z_msf = @tsm_membership_Gauss;
                otherwise
                    error( 'TSModel/initialize: unknown MSF <%s>', opts.msf)
            end
            
            obj.LS_type = opts.method;
            
            switch obj.Type
                case 'Static'
                    
                    % Scheduling var = u / forget ProductSpace
                    obj.z = obj.z_fct( obj );
                    obj.zmu = obj.z_msf( obj.z(:,1:obj.nu), obj.v, obj.m );
                    obj.z = [ obj.z(:,1:obj.nu), ones(obj.N_ident,1) ];
                    
                    if strcmp( opts.method, 'local' ) % local LS
                        obj.Theta = [];
                        for iv = 1 : obj.nv
                            Phi = bsxfun( @times, obj.z, obj.zmu(:,iv)  );
                            theta = Phi \ obj.y_ident;
                            obj.Theta = [obj.Theta, theta ];
                            phi =  transpose( reshape(  theta, obj.nu+1,1 ) );
                            obj.B( iv, : ) = phi( 1:obj.nu );
                            obj.C( iv, 1 ) = phi( obj.nu+1 );
                        end
                        
                    else % global LS
                        
                        Phi = bsxfun(@times, obj.z, obj.zmu(:,1) );
                        for iv = 2 : obj.nv
                            Phi = [ Phi, bsxfun(@times, obj.z, obj.zmu(:,iv) ) ];
                        end
                        obj.Theta = Phi \ obj.y_ident;
                        phi = transpose(  reshape( obj.Theta,obj.nu+1,obj.nv) );
                        obj.B = phi( :, 1:obj.nu );
                        obj.C = phi( :, obj.nu+1 );
                        
                    end
                    
                case { 'ARX', 'OE' }
                    
                    % Scheduling var z = [ y(lags_y) | u(lags_u) | 1 ]
                    obj.z = obj.z_fct( obj, obj.u_ident );
                    n = size( obj.z,1);
                    obj.zmu = obj.z_msf( obj.z, obj.v, obj.m );
                    obj.z = [ obj.z, ones(n,1) ];

                    %ToDo: #z_lags <> #x_lags???
                    if  strcmp( opts.method, 'global' ) % local LS
                        
                        Phi = bsxfun(@times, obj.z, obj.zmu(:,1) );
                        for iv = 2 : obj.nv
                            Phi = [ Phi, bsxfun(@times, obj.z, obj.zmu(:,iv) ) ];
                        end
                        obj.Theta = Phi \ obj.y_ident(obj.z_maxlag+1:end);
                        phi = transpose(  reshape( obj.Theta ,obj.nA+obj.nB+1,obj.nv) );
                        obj.A = phi( :, 1:obj.nA );
                        obj.B = phi( :, obj.nA+(1:obj.nB) );
                        obj.C = phi( :, obj.nA+obj.nB+1 );

                    elseif  strcmp( opts.method, 'local' )

                        for iv = 1 : obj.nv
                            obj.Phi = bsxfun( @times, obj.z, obj.zmu(:,iv)  );
                            phi =  transpose( reshape(  obj.Phi \ obj.y_ident(1:n), obj.nA+obj.nB+1,1 ) );
                            obj.A(iv,:) = phi( 1:obj.nA );
                            obj.B(iv,:) = phi( obj.nA+(1:obj.nB) );
                            obj.C(iv,:) = phi( obj.nA+obj.nB+1 );
                        end
                        
                    else
                        error('TSModel/initialize: unknow LS initi <%s>', obj.method)
                    end
            end
        end
        
        function [ A,B,C ] = getLM( obj )
            % Get local model parameters
            
            A = obj.A;
            B = obj.B;
            C = obj.C;
        end
        
        function obj = setLM( obj, A,B,C )
            % Set local model parameters A,B,C
            
            if ~isempty(A) && ~isequal( size(A), [obj.nv,length(obj.x_lag_y)] )
                error( 'TSModel/setLM: dim A <> nv x lag_y' )
            end
            obj.A = A;
            obj.nA = size(A,2);
            
            if ~isequal( size(B), [obj.nv,length(obj.x_lag_u)] )
                error( 'TSModel/setLM: dim B <> nu x nv x lag_u' )
            end
            
            obj.B = B;
            obj.nB = size(B,2);
            if ~isequal( size(C), [obj.nv,1] )
                error( 'TSModel/setLM: dim C <> nv x 1' )
            end
            obj.C = C;
        end
        
        function obj = set_msf_type( obj, type )
            % Set type of MSF: FCM or Gauss
            
            switch type
                case 'FCM'
                    obj.z_msf = @tsm_membership_FCM;
                case 'Gauss'
                    obj.z_msf = @tsm_membership_Gauss;
                otherwise
                    error( 'tsModel/set_msf_type: unknown type <%s>',type)
            end
        end
        
        function yp = predict( obj, u, y )
            % Predict the TS model at vectors u (OE) or u,y (ARX)
            
            % Static: len(u) = len(y)  = N
            % ARX: len(u) = len(y) = N
            % OE:  len(u) = N, y_0: len(y)>=max(z_lag_y,x_lag_y)
            
            % Check A/B/C not empty (model not yet initialized)
            if ~strcmp( obj.Type, 'Static') && isempty(obj.A)
                error( 'tsModel/predict: matrix A empty' )
            end
            if isempty(obj.B)
                error( 'tsModel/predict: matrix B empty' )
            end
            if isempty(obj.C)
                error( 'tsModel/predict: matrix C empty' )
            end
            
            switch obj.Type
                case 'Static'
                    yp = tsm_predict_Static( obj, u );
                case 'ARX'
                    yp = tsm_predict_ARX( obj, u, y );
                case 'OE'
                    yp = tsm_predict_OE( obj, u, y );
            end
        end
        
        function yp = simulate( obj, u, y )
            % Simulate the TS model with vectors u and y_0
            
            % Static: len(u) = len(y)
            % ARX: len(u)=N, len(y)=N
            % OE:  len(u)=N, y_0: len(y)=max(z_lag_y,x_lag_y)

            warning( 'TSModel/simulate: not yet implemented' )
            switch obj.Type
                case 'Static'
                    yp = tsm_predict_Static( obj, u );
                case 'ARX'
                    yp = tsm_predict_ARX( obj, u, y );
                case 'OE'
                    yp = tsm_predict_OE( obj, u, y );
            end
        end
        
        function obj = optimize( obj, method, varargin )
            % Optimize the TS model parameters
            
            p = inputParser;
            p.addRequired( 'method', @ischar )
            
            valOpt = @(x) isa(x,'optim.options.Lsqnonlin');
            p.addParameter('optimopts',[], valOpt)
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
            
            tic
            switch obj.Type
                
                case 'Static'
                    nvu = obj.nv*obj.nu;
                    switch upper(method(1))
                        case 'M' % optimize membership MF
                            % Limit v to range of u|y
                            obj.Opt_Type = 'MF';
                            lb = reshape( repmat(obj.Limits(1:obj.nu,1),1,obj.nv),1,nvu );
                            ub = reshape( repmat(obj.Limits(1:obj.nu,2),1,obj.nv),1,nvu );
                            x0 = [ reshape(obj.v,1,nvu) ];
                            [x_opt,resnorm,residual,exitflag,output] = ...
                                lsqnonlin( @(x) tsm_optimize_Static_MF( x,...
                                obj.u_ident, obj )-obj.y_ident,...
                                x0, lb,ub,optimopts);
                            obj.v = reshape( x_opt(1:nvu), obj.nv, obj.nu );
                        case 'L' % optimize local models LM
                            obj.Opt_Type = 'LM';
                            lb = []; % no constraints for B/c
                            ub = [];
                            x0 = [ reshape(obj.B,1,nvu),...
                                   reshape(obj.C,1,obj.nv) ];
                            [x_opt,resnorm,residual,exitflag,output] = ...
                                lsqnonlin( @(x) tsm_optimize_Static_LM( x,...
                                obj.u_ident, obj )-obj.y_ident,...
                                x0, lb,ub,optimopts);
                            n1 = 1; n2 = nvu;
                            obj.B = reshape( x_opt(n1:n2), obj.nv,obj.nu );
                            n1 = n2 + 1; n2 = n1+obj.nv-1;
                            obj.C = reshape( x_opt(n1:n2),obj.nv,1 );
                        case 'B' % optimize MF and LM
                            obj.Opt_Type = 'MF&LM';
                            lb = [ reshape( repmat(obj.Limits(1:obj.nu,1),1,obj.nv),1,nvu ),...
                                   -inf*ones(1,obj.nv*(obj.nu+1)) ];
                            ub = [ reshape( repmat(obj.Limits(1:obj.nu,2),1,obj.nv),1,nvu ),...
                                   +inf*ones(1,obj.nv*(obj.nu+1)) ];
                            x0 = [ reshape(obj.v,1,nvu),...
                                   reshape(obj.B,1,nvu),...
                                   reshape(obj.C,1,obj.nv) ];
                            [x_opt,resnorm,residual,exitflag,output] = ...
                                lsqnonlin( @(x) tsm_optimize_Static( x,...
                                obj.u_ident, obj )-obj.y_ident,...
                                x0, lb,ub,optimopts);
                            n2 = nvu;
                            obj.v = reshape( x_opt(1:n2), obj.nv, obj.nu );
                            n1 = n2 + 1; n2 = n1 + nvu - 1;
                            obj.B = reshape( x_opt(n1:n2), obj.nv,obj.nu );
                            n1 = n2 + 1; n2 = n1 + obj.nv-1;
                            obj.C = reshape( x_opt(n1:n2),obj.nv,1 );
                        otherwise
                            error('TSModel/optimize: method not M, L, or B' )
                    end
                    obj.OptimPR.optimopts = optimopts;
                    obj.OptimPR.exitflag = exitflag;
                    obj.OptimPR.resnorm = resnorm;
                    obj.OptimPR.residual = residual;
                    obj.OptimPR.output = output;
                    
                case 'ARX'
                    switch method
                        case 'M' % optimize membership MF
                            % Limit v to range of u|y
                            obj.Opt_Type = 'MF';
                            nu = obj.nB * obj.nu; % Lags beachten
                            lb = reshape( repmat(obj.Limits(1:obj.nu,1),1,obj.nv),...
                                1,obj.nv*obj.nu );
                            ub = reshape( repmat(obj.Limits(1:nu,2),1,obj.nv),1,obj.nv*obj.nu );
                            x0 = [ reshape(obj.v,1,obj.nv*obj.nu) ];
                            [x_opt,resnorm,residual,exitflag,output] = ...
                                lsqnonlin( @(x) tsm_optimize_ARX_MF( x,...
                                obj.u_ident, obj )-obj.y_ident,...
                                x0, lb,ub,optimopts);
                            obj.v = reshape( x_opt(1:obj.nv*obj.nu), obj.nv, obj.nu );
                        case 'L' % optimize local models LM
                            obj.Opt_Type = 'LM';
                            lb = [];
                            ub = [];
                            x0 = [ reshape(obj.B,1,obj.nv*obj.nu),...
                                reshape(obj.C,1,obj.nv) ];
                            [x_opt,resnorm,residual,exitflag,output] = ...
                                lsqnonlin( @(x) tsm_optimize_ARX_LM( x,...
                                obj.u_ident, obj )-obj.y_ident,...
                                x0, lb,ub,optimopts);
                            n1 = 1; n2 = obj.nv * obj.nu;
                            obj.B = reshape( x_opt(n1:n2), obj.nv,obj.nu );
                            n1 = n2 + 1; n2 = n1+obj.nv-1;
                            obj.C = reshape( x_opt(n1:n2),obj.nv,1 );
                        case 'B' % optimize MF and LM ( c | A | B C ]
                            obj.Opt_Type = 'MF&LM';
                            % lag_y | lag_u
                            % v: z_lag_u | zlag_y
                            lb = [];
                            ub = [];
                            for iu=1:obj.z_lag_nu
                                lb = [ lb, repmat( obj.Limits(iu,1),1,obj.nv*length(obj.z_lag_u{iu}) ) ];
                                ub = [ ub, repmat( obj.Limits(iu,2),1,obj.nv*length(obj.z_lag_u{iu}) ) ];
                            end
                            lb = [ lb, repmat( obj.Limits(obj.nu+1,1),1,obj.nv*length(obj.z_lag_y) ),...
                                      -inf*ones(1,obj.nv*(obj.nA+obj.nB+1)) ];
                            ub = [ ub, repmat( obj.Limits(obj.nu+1,2),1,obj.nv*length(obj.z_lag_y) ),...
                                       +inf*ones(1,obj.nv*(obj.nA+obj.nB+1)) ];
                            p0 = [ reshape(obj.v,1,obj.nv*(obj.nA+obj.nB)),...
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
                            
                            n2 = obj.nv * (obj.nA+obj.nB);
                            obj.v = reshape( p_opt(1:n2), obj.nv, (obj.nA+obj.nB) );
                            obj.Theta = transpose( p_opt(n2+1:end) );

                            n1 = n2 + 1; n2 = n1 + obj.nv*obj.nA-1;
                            obj.A = reshape( p_opt(n1:n2), obj.nv,obj.nA );
                            n1 = n2 + 1; n2 = n1 + obj.nv*obj.nB-1;
                            obj.B = reshape( p_opt(n1:n2), obj.nv,obj.nB );
                            n1 = n2 + 1; n2 = n1+obj.nv-1;
                            obj.C = reshape( p_opt(n1:n2),obj.nv,1 );
                        
                        otherwise
                            error('TSModel/optimize: not M, L, or B' )
                    end
                    obj.OptimPR.optimopts = optimopts;
                    obj.OptimPR.exitflag = exitflag;
                    obj.OptimPR.resnorm = resnorm;
                    obj.OptimPR.residual = residual;
                    obj.OptimPR.output = output;
                    
                case 'OE'
                       obj.Opt_Type = 'MF&LM';
                            % c = [lag_y,lag_u)
                            %ToDo: u x nu
                             lb = [ repmat( ...
                                    [ repmat(obj.Limits(2,1),1,obj.nA),...   % lag_y
                                      repmat(obj.Limits(1,1),1,obj.nB) ],...  % nu * lag_u
                                      1,obj.nv),... 
                                   -inf*ones(1,obj.nv*(obj.nA+obj.nB+1)) ];
                             ub = [ repmat( ...
                                    [ repmat(obj.Limits(2,2),1,obj.nA),...    % lag_y
                                      repmat(obj.Limits(1,2),1,obj.nB) ],...  % nu * lag_u
                                      1,obj.nv),... 
                                   +inf*ones(1,obj.nv*(obj.nA+obj.nB+1)) ];
                            
                            p0 = [ reshape(obj.v,1,obj.nv*(obj.nA+obj.nB)),...
                                   transpose(obj.Theta) ];
                            
                            [p_opt,resnorm,residual,exitflag,output] = ...
                                lsqnonlin( @(p) tsm_optimize_OE( p, obj ),...
                                p0, lb,ub,optimopts );
                            n2 = obj.nv * (obj.nA+obj.nB);
                            obj.v = reshape( p_opt(1:n2), obj.nv, (obj.nA+obj.nB) );
                            obj.Theta = transpose( p_opt(n2+1:end) );

                            n1 = n2 + 1; n2 = n1 + obj.nv*obj.nA-1;
                            obj.A = reshape( p_opt(n1:n2), obj.nv,obj.nA );
                            n1 = n2 + 1; n2 = n1 + obj.nv*obj.nB-1;
                            obj.B = reshape( p_opt(n1:n2), obj.nv,obj.nB );
                            n1 = n2 + 1; n2 = n1+obj.nv-1;
                            obj.C = reshape( p_opt(n1:n2),obj.nv,1 );
                        
            end
            obj.debug(2,sprintf('Duration of optimization: %g s\n', toc))

        end
        
        function obj = setName( obj, name )
            % Set name of TS model
            obj.Name = name;
        end
       
        function [A,B,C] = Regressor2ABC( obj, Theta )
            % Get local model matrices A,B,C from regressor matrix

            if nargin < 2
                Theta = obj.Theta;
            end
            n2 = obj.nv * obj.nA;
            A = reshape(Theta(1:n2), obj.nv, obj.nA);
            n1 = n2+1; n2 = obj.nv * obj.nB;
            B = reshape(Theta(n1:n2), obj.nv, obj.nB);
            n1 = n2+1; n2 = obj.nv * obj.nC;
            C = reshape(Theta(n1:n2), obj.nv, 1);
        end
        
        function disp( obj )
            % Dsiplay TS model settings and parameters
            
            fprintf( 'TS-Model: Type=%s\n', obj.Type )
            
            disp@tsm_Base( obj )

            fprintf( ' Structural parameters: nu = %d, ny = %d, nv = %d\n', obj.nu, obj.ny, obj.nv )
            
            if obj.N_ident > 0
                fprintf( ' Identification data: N=%d\n', obj.N_ident )
                if obj.ts_ident > 0
                    fprintf( ', ts=%g', obj.ts_ident )
                end
                if obj.C_ident
                    fprintf( ' (%s)', obj.C_ident )
                end
            end

            fprintf( ' Initial model estimation:\n')
            if ~isempty( obj.z_lag_u )
                fprintf( ' lags: ')
                for i=1:obj.nu
                    fprintf( 'u_%d:%s, ', i,mat2str( obj.z_lag_u{i} ) )
                end
            end
            
            if ~isempty( obj.z_lag_y )
                fprintf( 'y = %s\n', mat2str(obj.z_lag_y) )
            end
            
            if ~isempty( obj.z_lag_y )
                fprintf( '  Membership function type = %s\n', obj.z_Type )
            end
            switch obj.v_Type 
                case'FCM'
                    fprintf( '  Clustering: %s, nue=%g norm=%s', obj.v_Type, obj.nue, obj.v_Norm )
                case'Gauss'
                    fprintf( '  Clustering: %s, sigma=%g', obj.v_Type, obj.sigma )
            end
            if obj.ProductSpace
                    fprintf( ' in product space\n')
            else
                    fprintf( ' in input space\n')
            end
                
            fprintf( ' Estimation of local models:\n')
            if ~isempty( obj.x_lag_u )
                fprintf( ' lags:  ')
                for i=1:obj.nu
                    fprintf( '  u_%d:%s, ', i,mat2str( obj.x_lag_u{i} ) )
                end
            end
            if ~isempty( obj.z_lag_y )
                fprintf( '  y = %s\n', mat2str(obj.x_lag_y) )
            end
            if obj.LS_type
                fprintf( '  Initialization of local models: %s\n', obj.LS_type )
            end
            
            if obj.Opt_Type
                fprintf( '  Optimization of model parameters: %s\n', obj.Opt_Type )
            end
        end
        
        function h = plotCluster( obj, v, varargin )
            % Plot clustering as 2D-plots for all dimensions u and y
            p = inputParser;
            p.addRequired('v',@ismatrix)
            p.addParameter('figure',2,@isscalar)
            p.addParameter('title','Cluster',@ischar)
            p.addParameter('file','',@ischar)
            p.parse( v, varargin{:} )
            opts = p.Results;
            
            if nargin < 2
                v = obj.v;
            end
            n = size( v, 2 );
            
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
                        plot( obj.z(:,i1),obj.z(:,i2),'k.', 'MarkerSize',8 )
                        hold on
                        l{end+1} = 'data-points';
                    end
                    plot( v(:,i1), v(:,i2), 'rx', 'MarkerSize', 12 )
                    if ~isempty( obj.Limits )
                        axis( [obj.Limits(i1,:), obj.Limits(i2,:) ] )
                    end
                    axis square
                    grid on
                    box on
                    set( gca, 'FontSize', 14 )
                    if ~isempty( obj.Labels )
                        xlabel( obj.Labels(i1), 'FontSize', 14  )
                        ylabel( obj.Labels(i2), 'FontSize', 14  )
                    else
                        xlabel( sprintf('u_{%d}',i1), 'FontSize', 14  )
                        ylabel( sprintf('u_{%d}',i2), 'FontSize', 14  )
                    end
                    if is == 1 && ~isempty( opts.title )
                        title( opts.title, 'FontSize', 14 )
                        l{end+1} = 'cluster centers';
                        legend( l, 'location','SE', 'FontSize', 14  )
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
        
        function h = plotIdentData( obj, t,u,y,fig  )
            % Plot the identification data    
   
            
            if nargin < 4
                t = obj.t_ident;
                u = obj.u_ident;
                y = obj.y_ident;
                fig = 1;
            end

            %ToDo: u separat/zusammen
            h = figure( fig );
            clf

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
                    title( sprintf( 'Identification data N=%d data-points', obj.N_ident ) )
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
            % Plot the TS model: data and cluster centers
            
            plotIdentData( obj, obj.t_ident, obj.u_ident,  obj.y_ident );
            plotCluster( obj, obj.v );
            
        end
        
    end
    
end