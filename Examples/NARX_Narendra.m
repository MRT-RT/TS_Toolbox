%% TS-Toolkbox: ARX example (Narenda function)
%
% $Id$

%%
% Aproximate the Narendra function
%%
% 
% $$y_{k+1} = \frac{ y_k\cdot y_{k-1} \cdot ( y_k + 2.5 ) }{ 1 + y_k^2 + y_{k-1}^2 } +  u_k$$
%
% with  $y_1 = 0, y_2 = 0$
%
% as a NARX/TS with 3 local models
%%
% Source; K. Narendra and K. Parthasarathy. "Identification and control of dynamical systems using neural networks". IEEE Transactions on Neural Networks 1(1) (Mar. 1990), pp. 4-27.
%
% <https://www.sciencedirect.com/science/article/pii/0888613X9290014Q>

clear, close all

%%
% Path to TSModel files
addpath( '..', '../Functions' ) 

%% Model order
% number of local models
nc = 3;    
%%
% number of inputs
nu = 1;   
%%
% number of outputs
ny = 1;   
%%
% fuzziness parameter
nue = 1.2; 
%%
% sampling time
dt = 1e-2; 

%% Generate Narendra model data
%
% Use n data points
n = 1000;
rng(0);
[ u, y ] = Narendra_fct( n );

%%
% Generate time vector
t = dt * transpose( 0:size(u,1)-1 ); 


%% Plot of Narendra function

figure(1),clf
subplot(2,1,1)
plot(t,u)
grid
ylabel('u')
title(sprintf('Narendra function / u=steps / n=%d',n))

subplot(2,1,2)
plot(t,y)
grid
ylabel('y')
xlabel('t')

%% Create the TS model
ts = TSModel( 'ARX', nc, nu, 'Name','NARX', 'Comment','Narendra');
%
% Set lags: u(t-0), y(t-1),y(t-2)
ts.setLags( [0], [1,2] );
%%
% set tainging data u,y
ts.setData( u, y, 'SampleTime',dt, 'Labels', { 'u', 'y' } );
%%
% set data limis: u=[-2,2], y=[-5,10]
ts.setDataLimits( [-2,2 ; -5,10] );

%% Clustering in product-space (u,y)
ts.clustering( 'FCM', 'nue', nue, 'tries',1, 'seed', 0 )
c1 = getCluster(ts);

%% Initialize local models with LS (global/local)
ts.initialize( 'FBF', 'nue', nue, 'method','global'  );

%% Compute TS model for given data: yp
yp = ts.evaluate( u,y );

%%
% Plot residuals
plotResiduals( y, yp, 'figure', 2, 'title', 'Residuals Narendra/NARX' );

%% Plot training data y and predicted data yp
figure(3),clf

yyaxis left
plot(t,u,'b--')
ylabel('u')

yyaxis right
plot(t,y,'g-',t,yp,'r-')

grid on
ylabel('y')
xlabel('t')
title( sprintf('Narendra: ARX model ident c=%d \\nu=%g',ts.nc,ts.nue))
legend('u','y','y_{pred}')

%% Evaluate with test data ut,yt
[ut,yt] = Narendra_fct( n );
ypt = ts.evaluate( ut,yt );

%%
% Residual
plotResiduals( y, ypt, 'figure', 4, 'title', 'Residuals Narendra eval/NARX' );

figure(5),clf

yyaxis left
plot(t,ut,'b--')
ylabel('u')

yyaxis right
plot(t,yt,'g-',t,ypt,'r-')

grid on
ylabel('y')
xlabel('t')
title( sprintf('Narendra: ARX model eval c=%d \\nu=%g',ts.nc,ts.nue))
legend('u','y','y_{pred}')

%% Optimize Clusters c (MF) and/or local model A/B/C

ts.optimize( 'B' )
c2 = getCluster( ts );

%%
% Plot the residuals
ypo = ts.evaluate( u,y );
plotResiduals( y, ypo, 'figure', 6, 'title', 'Residuals Narendra/NARX opt' );

%% 
% Plot original function and optimized model
figure(7),clf

yyaxis left
plot(t,u,'b--')
ylabel('u')

yyaxis right
plot(t,y,'g-',t,ypo,'r-')

grid on
ylabel('y')
xlabel('t')

title( sprintf('Narendra: NARX model opt c=%d \\nu=%g',ts.nc,ts.nue))
legend('u','y','y_{pred}')
