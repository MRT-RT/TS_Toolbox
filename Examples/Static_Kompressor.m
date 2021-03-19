%% Takagi-Sugeno Model Identification Toolbox
%
% Static LiP model for a compressor model.
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
% Add path to TS model class
addpath( '../TSModel' )

%% Identification data 
% Input vector u and output vector y
load( 'Data/Kompressor.mat' );
u_val = [ KompDaten.x1, KompDaten.x2 ];
y_val = KompDaten.y;
dt = 1; % Implicit sampling time for static models


%% Structural parameters
nc = 6;       % number of clusters = local models
nu = 2;       % number of inputs
nue = 1.2;    % fuzziness parameter
MSF  = 'FCM'; % FCM type membership functions

%% Estimatiuon of TS model
ts = TSModel( 'Static', nc, nu, 'Name','LiP','Comment','Compressor');

ts.setData( u_val, y_val, 'SampleTime',dt, 'Labels', { 'p_R','eta', 'dm/dt' } );
ts.setDataLimits( [0.2,1 ; 1,2.2; 0,200] );

%% 
% Clustering in product-space: $[u,y]$
ts.clustering( 'FCM', 'nue', nue, 'tries', 10, 'seed', 0 );
v1 = getCluster( ts );

%% 
% Initialization of local LS models with global Least-squares
ts.initialize( 'FCM', 'nue', nu, 'method','global'  );

%% 
% Predict the TS model for observed input data
y_pred = ts.predict( u_val );
plotResiduals( y_val, y_pred, 'figure', 4, 'title', 'Compressor: correlation' );

%% 
% Optimize TS model: optimize only cluster centers $v$ 
ts.optimize( 'M' );
v2 = getCluster( ts );

%%
% Plot the correlation
y_pred_opt = ts.predict( u_val );
plotResiduals( y_val, y_pred_opt, 'figure', 5, 'title', 'Compressor opt: correlation' );

% Plot the correlation histogram
hv = plotResidualHist( y_val, y_pred_opt, 'nbins',21, ...
        'figure', 6, 'title', 'Compressor opt: correlation histogram','nbins',21 );
hv.WindowState = 'maximized';

