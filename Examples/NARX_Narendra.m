%% Takagi-Sugeno Model Identification Toolbox
%
% NARX LiP model example for the Narendra function
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

% $Id$

%%
% Aproximate the Narendra function
%%
% 
% $$ y_{k+1} = \frac{ y_k\cdot y_{k-1} \cdot ( y_k + 2.5 ) }{ 1 + y_k^2 + y_{k-1}^2 } +  u_k$$
%
% with  inital states $y_1 = 0, y_2 = 0$
%
% as a NARX/TS with $n_v=3$ local models
%%
% Source; K. Narendra and K. Parthasarathy. "Identification and control of dynamical systems using neural networks". IEEE Transactions on Neural Networks 1(1) (Mar. 1990), pp. 4-27.
%
% <https://www.sciencedirect.com/science/article/pii/0888613X9290014Q>

%% Model order
nc = 3; % number of local models    
nu = 1; % number of inputs   
ny = 1; % number of outputs  

%%
% Fuzziness parameter
nue = 1.2; 


%% Generate data from Narendra function:
%
% Use $N$ data points
rng(0); % Initalize random number generator
N = 1000;
[ u, y ] = Narendra_fct( n );

%%
% Sampling time
dt = 1e-2; 
%%
% Create time vector: $t$
t = dt * transpose( 0:size(u,1)-1 ); 


%% Plot of the Narendra function
figure(1),clf

subplot(2,1,1)
plot(t,u)
grid on
ylabel('u')
title(sprintf('Narendra function / input: steps / N=%d',n))

subplot(2,1,2)
plot(t,y)
grid on
ylabel('y')
xlabel('t')

%% Create the TS model
ts = TSModel( 'ARX', nc, nu, 'Name','NARX Narendra', 'Comment','Narendra function');
%
% Set lags: u(t-0), y(t-1),y(t-2)
ts.setLags( {1}, [1,2] );
%%
% Set identification data $u$,$y$
ts.setData( u, y, 'SampleTime',dt, 'Labels', { 'u', 'y' } );
%%
% Set data limis: u=[-2,2], y=[-5,10]
ts.setDataLimits( [-2,2 ; -5,10] );

%% Clustering in product-space: $z=[u,y]$ with FCM membership functions:
ts.clustering( 'FCM', 'nue',nue, 'tries',1, 'seed',0 )
%% Get the found cluster centers
v1 = getCluster(ts);

%% Initialize the local models with global LS:
ts.initialize( 'FCM', 'nue', nue, 'method','global'  );

%% Predict the  TS model ouput for given data: $y_{pred}$
y_pred = ts.predict( u, y );

%%
% Plot the residuals
plotResiduals( y, y_pred, 'figure', 2, 'title', 'Correlation Narendra/NARX' );

%% Plot training data y and predicted data y_pred
figure(3),clf

yyaxis left
plot(t,u,'b--')
ylabel('u')

yyaxis right
plot(t,y,'g-',t,y_pred,'r-')

grid on
ylabel('y')
xlabel('t')
title( sprintf('Narendra: ARX model ident c=%d \\nu=%g',ts.nv,ts.nue))
legend('u','y_{obsv}','y_{pred}')

%% Predict output from new test data $u_t$, $y_t$
[ut,yt] = Narendra_fct( n );
yp_t = ts.predict( ut,yt );

%%
% Residual
plotResiduals( y, yp_t, 'figure', 4, 'title', 'Residuals Narendra eval/NARX' );

figure(5),clf

yyaxis left
plot(t,ut,'b--')
ylabel('u')

yyaxis right
plot(t,yt,'g-',t,yp_t,'r-')

grid on
ylabel('y')
xlabel('t')
title( sprintf('Narendra: NARX model n_v=%d \\nu=%g',ts.nv,ts.nue))
legend('u','y_{obsv}','y_{pred}')

%% Optimize Clusters c (MF) and/or local model A/B/C

ts.optimize( 'B' )
v2 = getCluster( ts );

%%
% Plot the residuals
yp_o = ts.predict( u,y );
plotResiduals( y, yp_o, 'figure', 6, 'title', 'Residuals Narendra/NARX opt' );

%% 
% Plot original function and optimized model
figure(7),clf

yyaxis left
plot(t,u,'b--')
ylabel('u')

yyaxis right
plot(t,y,'g-',t,yp_o,'r-')

grid on
ylabel('y')
xlabel('t')

title( sprintf('Narendra: NARX model opt n_v=%d \\nu=%g',ts.nv,ts.nue))
legend('u','y_{obsv}','y_{pred}')

