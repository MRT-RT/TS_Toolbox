%% TS-Toolbox: NARX example (throttle valve)

% $Id$

clear, close all

addpath( '..' )  % Path to TSModel files

%% Model order
nc = 3;    % number of clusters = local models
nu = 1;    % number of inputs
ny = 1;    % number of outputs

nue = 1.1; % fuzziness parameter
 
% Input vector u and output vector y
load( 'Data/Throttle1.mat' );

dt = 1e-2; % Sampling time
t = dt * transpose( 0:size(u,1)-1 );

figure(1),clf
subplot(2,1,1)
plot(t,u)
grid
ylabel('u')
title( 'Throttle data' )
subplot(2,1,2)
plot(t,y)
grid
ylabel('y')
xlabel('t')

%% Create TS model
ts = TSModel( 'ARX', nc, nu, 'Name','ARX', 'Comment','Throtte valve')
ts.setSchedulingLags( [0], [1,2,3] );
ts.setRegressorLags( [0], [1,2,3] );

ts.setData( u, y, 'SampleTime',dt, 'Labels', { 'volt', 'angle' } );
ts.setDataLimits( [-15,55 ; 0,5] );

%% Clustering in product-space (u,y)
ts.clustering( 'FCM', 'nue', nue, 'tries',10, 'seed', 0 )
c1 = getCluster(ts);

%% Initialize local LS models with LS (global/local)
ts.initialize( 'FBF', 'nue', nue, 'method','global'  );

%% Compute TS model for given data
yp = ts.evaluate( u,y );
plotResiduals( y, yp, 'figure', 2, 'title', 'Residuals Throttle NARX' );

return

%% Optimize Clusters c (MF) and/or local model A/B/C
ts.optimize( 'B' )
c2 = getCluster( ts );
ypo = ts.evaluate( u,y );
plotResiduals( y, ypo, 'figure', 3, 'title', 'Residuals Throttle NARX opt' );

