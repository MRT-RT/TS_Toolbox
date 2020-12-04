%% TS-Toolox: LiP example (compressor data)

clear, close all

addpath( '../TSModel' ) % Path to TSModel files

%% Modell order
nc = 6;    % number of clusters = local models
nu = 2;    % number of inputs
nue = 1.2; % fuzziness parameter
 
% Input vector u and output vector y
load( 'Data/Kompressor.mat' );
u = [ KompDaten.x1, KompDaten.x2 ];
y = KompDaten.y;

dt = 1; % Implicit sampling time for static models

%% Create TS model
ts = TSModel( 'LiP', nc, nu, 'Name','LS','Comment','Kompressor')

ts.setData( u, y, 'SampleTime',dt, 'Labels', { 'p_R','eta', 'dm/dt' } );
ts.setDataLimits( [0.2,1 ; 1,2.2; 0,200] );

%% Clustering in product-space (u,y)
ts.clustering( 'FCM', 'nue', nue, 'tries',10, 'seed', 0 )
c1 = getCluster(ts);

%% Initialize local LS models with LS (global/local)
ts.initialize( 'FBF', 'nue', nu, 'method','global'  );

%% Compute TS model for given data
yp = ts.evaluate( u );
plotResiduals( y, yp, 'figure', 4, 'title', 'Residuals Compressor' );

%% Optimize Clusters c (MF) and/or local model A/B/C
ts.optimize( 'B' )
c2 = getCluster( ts );
ypo = ts.evaluate( u );
plotResiduals( y, ypo, 'figure', 5, 'title', 'Residuals Compressor opt' );

