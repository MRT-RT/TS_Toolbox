%% Takagi-Sugeno Model Identification Toolbox
%
% Example of a NOE LiP TS model for the Narendra function
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
% $Id$%

%%
% Path to TSModel files
addpath( '..' ) 

%% Model order
nc = 3;    % number of clusters = local models
nu = 1;    % number of inputs
ny = 1;    % number of outputs
nue = 1.2; % fuzziness parameter
%%
% Scheduling lags $z$ = regressor lags $x=[y(k-1),y(k-2),u(k)]$ 
z_lag_u = {0};
z_lag_y = [2];
x_lag_u = {0};
x_lag_y = [2];

%% Compute identification data
%
% Create input $u$ as steps with width $l=1,\ldots,20]$ for $N=1000$ time steps 
% and compute the output $y$ from the Narendra function 
N = 1000;
rng(0);
[u,y] = Narendra_fct( N );
%%
% Sampling time
dt = 1e-2; 
%%
% Creat tiume vector $t$
t = dt * transpose( 0:size(u,1)-1 );

%% Plot of the identification data
h=figure(1);clf

subplot(2,1,1)
plot(t,u)
grid
ylabel('u')
subplot(2,1,2)
plot(t,y)
grid
ylabel('y')
xlabel('t')

%% Creation of  TS model
ts = TSModel( 'OE', nc, nu, 'Name','OE Narendra', 'Comment','Narendra function');
ts.setSchedulingLags( z_lag_u, z_lag_y );
ts.setRegressorLags( x_lag_u, x_lag_y );
%%
% Set the identification data
ts.setData( u, y, 'SampleTime',dt, 'Labels', { 'u', 'y' } );
ts.setDataLimits( [-2,2 ; -5,10] );

%% Clustering 
% Clustering in product-space $z=[u,y]$ with FCM membership functions and
% $\nu=1.2$ with 3 multi-start tries and ficed initialized random number
% generator
ts.clustering( 'FCM', 'nue', nue, 'tries',3, 'seed', 0 )
%%
% Get the cluster centers of the inital model
v1 = getCluster(ts)

%% Initialization of local models 
% with global Least-Squares, FCM membership funtions and $\nu=1.2$
ts.initialize( 'FCM', 'nue', nue, 'method','global'  );

%% Predicted NOE TS model ouput
yp = ts.predict( u,y );
hr=plotResiduals( y, yp, 'figure', 2, 'title', 'Narendra NOE: correlation' );
hr.WindowState = 'maximized';

%%
% Plot of the observed vs. predicted outputs
h=figure(3);clf

plot(t,u,'k-',t,y,'g-',t,yp,'r-')
grid on
xlabel( 'time t' )
ylabel( 'output y' )
title('Narendra NOE: observed vs. predicted outputs')
legend('u','y_{obsv}','y_{pred}')
h.WindowState = 'maximized';

%% Prediction on validation data
[ut,yt] = Narendra_fct( N );
ypt = ts.predict( ut,yt );
%%
% Plot the correlation
hr=plotResiduals( y, ypt, 'figure', 4, 'title', 'Narendra NOE: corrleation on test data' );
h.WindowState = 'maximized';
%%
% Plot of observed vs. predicted validation data
h=figure(5);clf

yyaxis left
plot(t,ut,'b--')
ylabel( 'u' ) 
yyaxis right

plot(t,yt,'g-',t,ypt,'r-')
ylabel( 'y' ) 
xlabel( 't' )

grid on
title('Narendra NOE: validation data')
legend('u','y','y_{pred}')
h.WindowState = 'maximized';

%% Optimization of the TS model parameters
% Optimize the cluster centers $v$ (MF) and the local model parameters $A_i, B_i, c_i$
ts.optimize( 'B' )
%%
% Get the cluster centers of the optimized NOE TS model
v2 = getCluster( ts )
%%
% Plot the correlation on validation data
ypo = ts.predict( u,y );
hr=plotResiduals( y, ypo, 'figure', 6, 'title', 'Narendra NOE: correlation opt' );
hr.WindowState = 'maximized';
%%
% Plot the observed vs. the predicted validation data
h=figure(7);clf
yyaxis left
plot(t,ut)
ylabel( 'input u' )
yyaxis right
plot(t,yt,'g-',t,ypt,'r-')
grid on
ylabel( 'output y' )
xlabel( 'time t' )
title('Narendra NOE optimized: predicted validation output')
legend('u','y_{obsv}','y_{pred}')
h.WindowState = 'maximized';

