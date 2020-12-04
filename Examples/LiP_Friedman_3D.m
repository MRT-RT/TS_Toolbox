%% TS-Toolbox: LiP example (Friedman function 3-dim)

% $Id$

clear, close all

addpath( '../TSModel' ) % Path to TSModel files

%% Modell order
nc = 3;    % number of clusters = local models
nu = 3;    % number of inputs
nue = 1.2; % fuzziness parameter
 
% Input vector u and output vector y
n = 500;
u = rand(n,nu);
y = Friedman_fct( u, nu );

dt = 1; % Implicit sampling time for static models

%% Create TS model
ts = TSModel( 'LiP', nc, nu, 'Name','Friedman3','Comment','Friedman function 3-dimensional')

ts.setData( u, y, 'SampleTime',dt, 'Labels', { 'u_1','u_2','u_3', 'y' } );
%ts.setDataLimits( [0,1 ; 0,1 ; 0,10] );

%% Clustering in product-space (u,y)
ts.clustering( 'FCM', 'nue', nue, 'tries',10, 'seed', 0 )
c1 = getCluster(ts);
ts.plotCluster( c1,'figure', 3, 'file','Friedman-c1', 'title', 'initial clusters' )

%% Initialize local LS models with LS (global/local)
ts.initialize( 'FBF', 'nue', nu, 'method','local'  );

%% Compute TS model for given data
yp = ts.evaluate( u );
plotResiduals( y, yp, 'figure', 4, 'title', 'Residuals Friedman-3' );

%% Optimize Clusters c (MF) and/or local model A/B/C
ts.optimize( 'B' )
c2 = getCluster( ts );
yp = ts.evaluate( u );
plotResiduals( y, yp, 'figure', 5, 'title', 'Residuals Friedman-3' );

%% Save the TS model
save( 'TSM_Friedman.mat', 'ts' )

%% Evaluate TS model on regular grid
ngrid = 50;
[u1g,u2g,u3g] = meshgrid( linspace(0,1,ngrid),linspace(0,1,ngrid),linspace(0,1,ngrid) );
ugrid = [ reshape( u1g,ngrid*ngrid*ngrid,1), reshape(u2g,ngrid*ngrid*ngrid,1), reshape(u3g,ngrid*ngrid*ngrid,1) ];

ygrid = Friedman_fct( ugrid, nu );
ypgrid = ts.evaluate( ugrid );

plotResiduals( ygrid, ypgrid, 'figure', 10, 'title', 'Residuals Friedman-3 grid data' );

%% 3D plot of function
figure(11),clf

i=find(ugrid(:,3) == u3g(1));
subplot(2,2,1)
plot3( ugrid(i,1),ugrid(i,2),ygrid(i),'.')
grid on
xlabel( 'u_1' ),ylabel( 'u_2' ),zlabel( 'y' )

i=find(ugrid(:,2) == u2g(1));
subplot(2,2,2)
plot3( ugrid(i,1),ugrid(i,3),ypgrid(i),'.')
grid on
xlabel( 'u_1' ),ylabel( 'u_3' ),zlabel( 'y_p' )

i=find(ugrid(:,1) == u1g(1));
subplot(2,2,3)
plot3( ugrid(i,2),ugrid(i,3),ygrid(i),'.')
grid on
xlabel( 'u_2' ),ylabel( 'u_3' ),zlabel( 'y_p' )

