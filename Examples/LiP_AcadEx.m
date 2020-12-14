%% TS-Toolbox: LiP example (Academic Example)

% $Id$

%%
% Approximate an academic example
%%
% 
% $$y = \sum_{i=1}^3 a_{i,}1\cdot u_1 + a_{i,2}\cdot u_2 + a_{i,0}$$
% 
% with given local models $a_i$

clear, close all
addpath( '../TSModel' ) % Path to TSModel files

%% Model order
nc = 3;    % number of clusters = local models
nu = 2;    % number of inputs
nue = 1.2; % fuzziness parameter

% Input vector u and output vector y
load( 'Data/AcadEx.mat' )
dt = -1; % Implicit sampling time for static models

%% Data of the true model (not really needed)
%
% Initial = true clusters
c0 = [ 0.5, 0.5
       0.5, 1.5
       1.5, 1.0 ];
%%
% Initial local models = true models
a0 = [ -4,  4, -2
        4, -2, -4
        2,  1,  1 ];

%% Create true model as TSModel: ts0
ts0 = TSModel( 'LiP', nc, nu, 'Name','Academic', 'Comment','true academic example');
ts0.setCluster( c0 );
ts0.setFuziness( nue );
ts0.set_msf_type( 'FBF' );

%%
% Set the parameters a of local models (A,B,C)
ts0.setLM( [], a0(:,1:nu), a0(:,end) );

%% Create TS model
ts1 = TSModel( 'LiP', nc, nu, 'Name','AcadEx','Comment','Academic example')

%%
% Set training u,y data and limits of u,y
ts1.setData( u, y, 'SampleTime',dt, 'Labels', { 'u_1','u_2', 'y' }, 'Comment', 'from function' );
ts1.setDataLimits( [0,2 ;0,2; -7,7] ); % [u1,u2], [y]

%% 
% Set the inital clusters c0
ts1.setCluster( c0 );
ts1.plotCluster( c0,'figure', 1, 'file','AcadEx-c0', 'title', 'AcadEx: initial clusters c_0' )

%% Clustering in product-space (u,y) with KMeans algorithm
% with 10 tries and random generator seed 0
ts1.clustering( 'KMeans', 'nue', nu, 'tries',10, 'seed', 0 )
c11 = getCluster( ts1 );

%% Initialize TS as FBF and LiP with global LS 
ts1.initialize( 'FBF', 'nue', nue, 'method', 'G'  );
% Show the settings of the TS object
disp( ts1 )

%% Compute TS model for given data
yp1 = ts1.evaluate( u );
% and show the residuals
plotResiduals( y, yp1, 'figure', 2, 'title', 'AcadEx: residuals' );

%% Optimize Clusters c (MF) and/or local models A/B/C or (Both)
ts1.optimize( 'B' )
c12 = getCluster( ts1 );
yp1 = ts1.evaluate( u );
plotResiduals( y, yp1, 'figure', 3, 'title', 'AcadEx: residuals of opt. model' );

%% Create copy of the TS model
ts2 = copy( ts1 );
%%
% Use Gustaffson-Kessel (GK) algorithm for clustering
ts2.clustering( 'GK', 'nue', nue, 'tolerance',1e-4 )
ts2.initialize( 'FBF', 'nue', nue  );
yp2 = ts2.evaluate( u );

%% Optimize only memberships
ts2.optimize( 'M' );
e2 = ErrorCriteria(y,yp2);
% Get the local model parameters (y=B2*u+C2)
[~,B2,C2] = getLM( ts2 )

%% Create another copy of the TS model
ts3 = copy( ts2 );

%%
% Optimize only local models
ts3.optimize( 'L' );
yp3 = ts3.evaluate( u );
e3 = ErrorCriteria( y, yp3 );

%% Evaluate TS model on regular grid 50x50
ngrid = 50;
[u1g,u2g] = meshgrid( linspace(0,2,ngrid),linspace(0,2,ngrid) );
ugrid = [ reshape( u1g,ngrid*ngrid,1), reshape(u2g,ngrid*ngrid,1) ];

% true values
ygrid = ts0.evaluate( ugrid );
% predicted values
ypgrid= ts1.evaluate( ugrid );

plotResiduals( ygrid, ypgrid, 'figure', 4, 'title', 'Residuals AcadEx grid data' );

%% Plot 2x2 true vs. predicted model
figure(5),clf

subplot(2,2,1)
plot3( ugrid(:,1),ugrid(:,2),ygrid,'.')
grid on
xlabel( 'u_1' ),ylabel( 'u_2' ),zlabel( 'y' )
title('true model')

subplot(2,2,3)
plot3( ugrid(:,1),ugrid(:,2),ypgrid,'.')
grid on
xlabel( 'u_1' ),ylabel( 'u_2' ),zlabel( 'y_p' )
title('predicted model')

subplot(2,2,2)
plot3( ugrid(:,1),ugrid(:,2),ypgrid-ygrid,'.')
grid on
xlabel( 'u_1' ),ylabel( 'u_2' ),zlabel( 'y-y_p' )
title('difference true - prediction model' )

subplot(2,2,4)
plot3( ugrid(:,1),ugrid(:,2),ygrid,'k.', ugrid(:,1),ugrid(:,2),ypgrid,'r.')
grid on
xlabel( 'u_1' ),ylabel( 'u_2' ),zlabel( 'y,y_p' )
title('model & prediction')
legend('model','prediction','location','nw')


%% Plot contour of true model
yg = reshape(ygrid,ngrid,ngrid);
ypg = reshape(ypgrid,ngrid,ngrid);

figure(6),clf
hold on
contour(u1g,u2g,ypg,50)
%%
% plot centers of local models
plot3( ts0.c(:,1),ts0.c(:,2),zeros(nc,1),'rx' )
colorbar
