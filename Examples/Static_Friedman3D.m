%% Takagi-Sugeno Model Identification Toolbox
%
% Static LiP model for the 3-dimensional Friedman test function.
%
% V1.0
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

%% Identification data 
% Use the 3-dimensional Friedman function: 
nu = 3;
%%
% $$ y = 10\cdot\sin( \pi\cdot u_1 \cdot u_2 + 20\cdot(u_3-0.5)^2 $$
%%
% Choose the fuzziness parameter $\nu = 1.2$
nue = 1.2;
%%
% Choose the input matrix $u$ as random data with $N$ data-points: $u_{1,2}\in[0,1]$
N = 500;
u = rand( N, nu );
%%
% Compute the output vector $y$:
y = Friedman_fct( u, nu );

%% Structural parameters
% Number of clusters $n_v$ = number of local models
nv = 5;    
%%
% Membership function Type: FCM
MSF = 'FCM';

%% Estimation of the static LiP TS model
%
addpath( '../TSModel' );  % Path to TSModel class
ts = TSModel( 'Static', nv, nu, 'comment', 'Friedman 3D' );
%%
% Set the identification data: $u$, $y$
ts.setData( u, y );

%%
% Clustering:
%%
% * FCM: fuzziness parameter $\nu=1.2$ with Euclidian norm (default)
% * clustering in product-space
% * Multi-Start: 5 tries
ts.clustering( MSF, 'nue', nue, 'productspace', true, 'tries', 5 );
            
%%
% Initialisation of local models: global least squares estimation
ts.initialize( MSF, 'nue', nue, 'method', 'global' );

%%
% Optimization of both: membership and local model parameters
ts.optimize( 'Both' );

%%
% Show the resultiung TS model parameters:
disp( ts )

%%
% Plot the cluster centers: $v$
v = getCluster( ts )
ts.plotCluster( v, 'figure',1);

%%
% Predict the TS model output: $y_{pred}$
y_pred = ts.predict( u, y );

%%
% Plot the correlation
plotResiduals( y, y_pred, 'figure', 2, 'title', 'Friedman-3D: correleation' );
set(gcf,'WindowState','maximized');
%%
% Plot the residual histogram:
plotResidualHist( y, y_pred, 'figure', 3, 'nbins', 21, ...
    'title', 'Friedman-3D: residual histogram' );
set(gcf,'WindowState','maximized');

%%
% Plot the rule activation and input/output data:
plotRuleActivation( u,y_pred, ts, 'figure', 4 );
set(gcf,'WindowState','maximized');

%% Validation of the TS model
%
% Choose  another $N$ random inputs $[u_1,u_2,u_3]$
u_val = rand( N, nu );
y_obsv = Friedman_fct( u_val, nu );
%%
% Compute the output vector: $y_{pred}$
y_val_pred = ts.predict( u_val );
%%
% Plot the outputs
figure(5);clf
plot( 1:N, y_obsv, 'k-',1:N, y_val_pred, 'r--' )
grid on
xlabel('k')
ylabel('y')
title( 'Friedman-3D: Validation' )
legend( 'y_{obsv}','y_{pred}' )
set(gcf,'WindowState','maximized');
%%
% Plot the correlation 
plotResiduals( y_obsv, y_val_pred, 'figure', 6, ...
    'title', 'Friedman-3D: validation/correlation' );
set(gcf,'WindowState','maximized');
%%
% Plot the correlation histogram
plotResidualHist( y_obsv, y_val_pred, 'figure', 7, 'nbins', 31, ...
    'title', 'Friedman-3D: validation/correlation histogram' );
set(gcf,'WindowState','maximized');
