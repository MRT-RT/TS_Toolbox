%% Takagi-Sugeno Model Identification Toolbox
%
% Example of a throttle valve as a nonlinear ARX model LiP TS model (NARX)
%
%
% Axel Dürrbaum (<axel.duerrbaum@mrt.uni-kassel.de>)
%
% Department of Measurement and Control (MRT)
%
% Institute for System Analytics and Control (ISAC) 
%
% University of Kassel, Germany 
% (<http://www.uni-kassel.de/go/mrt>) 
%
% $Id$

%%
% Estimation of a NARX LiP TS model from the measurement data from a throttle
% valve
%%
% $$ y_{k+1} = \sum_{i=1}^{n_v} \mu_i(z) \cdot \left( \sum A_i\cdot x_y +
% B_i\cdot x_u + c_i \right) $$
% with
%%
% 
% * scheduling vector $z=[y_{k-1},y_{k-2},y_{k-3},u_k]$
% * FCM clustering in product-space with 3 tries
% * regression vector $x=[x_y,x_u]$, $x_y=[y_{k-1},y_{k-2},y_{k-3}]$ and $x_u=[u_k]$
% * initalization of the local models with a global Least-Squares
% * optimization of the parametes of the cluster centers $v$ and the local
% models $A,B,c$

%%
% Path to TSModel class
addpath( '../TSModel' )  

%% Structural parameters
nv = 3;    % number of clusters = local models
nu = 1;    % number of inputs
ny = 1;    % number of outputs
nue = 1.1; % fuzziness parameter
MSF = 'FCM';

%% Identification data 
%
% Load input vector $u$ and output vector $y$ from file:
load( 'Data/Throttle1.mat' );
%%
% Compute the time vector
dt = 1e-2; % Sampling time
t = dt * transpose( 0:size(u,1)-1 );

%% Creation of the TS model
ts = TSModel( 'ARX', nv, nu, 'Name','NARX', 'Comment','Throtte valve');
%%
% Set the identification data
ts.setData( u, y, 'SampleTime',dt, 'Labels', { 'volt', 'angle' } );
ts.setDataLimits( [-20,45 ; 0,5] );
%%
% Plot the identification data
hp=plotIdentData( ts );
hp.WindowState = 'maximized';

%%
% Set the scheduling lags:  u(k), y(k-1), y(k-2),  y(k-3)
ts.setSchedulingLags( [0], [1,2,3] );
%%
% Set the regressor lags:  u(k), y(k-1), y(k-2),  y(k-3)
ts.setRegressorLags( [0], [1,2,3] );

%% Clustering 
%
% Clustering is done in product-space (u,y) with FCM and $\nu=1.1$ for a
% multi-start of 10 tries and random number fixed initialized (Seed)
ts.clustering( 'FCM', 'nue', nue, 'tries',10, 'seed', 0 );
%%
% Cluster centers $v_1$: (columns $z=[y_{k-1},y_{k-2},y_{k-3},u_k]$, rows=local models)
v1 = getCluster(ts)

%% Initialization of the local models 
%
% The local models are intialized with global Least-Squares and FCM
% membership functions ($\nu=1.1$)
ts.initialize( MSF, 'nue', nue, 'method','global'  );

%% Prediction of the initial TS model 
%
% Prediction of the NARX TS model with the identificaton input: $u$
y_pred = ts.predict( u,y );
hi=plotResiduals( y, y_pred, 'figure', 2, 'title', 'Residuals Throttle NARX' );
hi.WindowState = 'maximized';

%% Optimization of the TS model
%
% by optimizing clusters centers $v$  (MF) and local model parameters $A/B/c$
ts.optimize( 'B' );
%%
% Cluster centers $v_1$: (columns: $z=[y_{k-1},y_{k-2},y_{k-3},u_k]$, rows=local models)
v2 = getCluster( ts )
%%
% Local model matrices $A_i\rightarrow x_y,B_i\rightarrow x_u, c_i$
[A,B,c] = getLM( ts )

%% Prediction of the optimized TS model
%
y_pred_opt = ts.predict( u,y );
ho=plotResiduals( y, y_pred_opt, 'figure', 3, 'title', 'Residuals Throttle NARX opt' );
ho.WindowState = 'maximized';
