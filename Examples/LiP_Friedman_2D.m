%% TS-Toolbox: LiP example (Friedman function 2-dim)

% $Id%

clear, close all

addpath( '../TSModel' ) % Path to TSModel files

%% Modell order
nc = 3;    % number of clusters = local models
nu = 2;    % number of inputs
nue = 1.2; % fuzziness parameter

% Input vector u and output vector y
n = 500;
u = rand(n,nu);
y = Friedman_fct( u, nu );

dt = 1; % Implicit sampling time for static models

%% Create TS model
ts = TSModel( 'LiP', nc, nu, 'Name','Friedman2', 'Comment','Friedman function 2-dimensional')

ts.setData( u, y, 'SampleTime',dt, 'Labels', { 'u_1','u_2', 'y' } );
%ts.setDataLimits( [0,1 ; 0,1 ; 0,10] );

%% Clustering in product-space (u,y)
ts.clustering( 'FCM', 'nue', nue, 'tries',10, 'seed', 0 )
c1 = getCluster( ts );
ts.plotCluster( c1,'figure', 3, 'file','Friedman-c1', 'title', 'initial clusters' )

%% Initialize local LS models with LS (global/local)
ts.initialize( 'FBF', 'nue', nu, 'method','local'  );

%% Compute TS model for given data
yp = ts.evaluate( u );
plotResiduals( y, yp, 'figure', 4, 'title', 'Residuals Friedman-2' );

%% Optimize Clusters c (M) or local models  A/B/C (L) or both (B)
ts.optimize( 'B' )
c2 = getCluster( ts );
yp = ts.evaluate( u );
plotResiduals( y, yp, 'figure', 5, 'title', 'Residuals Friedman-2' );

%% Evaluate TS model on regular grid
ngrid = 50;
[u1g,u2g] = meshgrid( linspace(0,1,ngrid),linspace(0,1,ngrid) );
ugrid = [ reshape( u1g,ngrid*ngrid,1), reshape(u2g,ngrid*ngrid,1) ];

ygrid = Friedman_fct( ugrid, 2);
ypgrid = ts.evaluate( ugrid );

plotResiduals( ygrid, ypgrid, 'figure', 10, 'title', 'Residuals Friedman-2 grid data' );

%% 3D plot of function
figure(11),clf
subplot(2,2,1)
plot3( ugrid(:,1),ugrid(:,2),ygrid,'.')
grid on
xlabel( 'u_1' ),ylabel( 'u_2' ),zlabel( 'y' )
subplot(2,2,3)
plot3( ugrid(:,1),ugrid(:,2),ypgrid,'.')
grid on
xlabel( 'u_1' ),ylabel( 'u_2' ),zlabel( 'y_p' )
subplot(2,2,2)
plot3( ugrid(:,1),ugrid(:,2),ypgrid-ygrid,'.')
grid on
xlabel( 'u_1' ),ylabel( 'u_2' ),zlabel( 'y-y_p' )

