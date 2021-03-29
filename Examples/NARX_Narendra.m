%% Takagi-Sugeno Model Identification Toolbox
%
% Example of a NARX TS model for a Narendra function
%
% V 1.0
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
% Example of the identification of a NARX MISO TS model for 
% given multiple inputs $u$ and single output $y$.
%% 
% Determine the NARX TS model
%
% $$ \hat{y}_{k+1} = \sum_{i=1}^{n_v} \phi_i(z) \cdot \left( \sum_{l=1}^{l_y} A_{i}\cdot y_{k-l} + \sum_{j=1}^{n_u}\sum_{l=0}^{l_u} B_{i,j}\cdot u_{k-l} + c_{i} \right)$$
%%
% * for given $u_j, j=1,\ldots,n_u$ of $n_u$ input vectors and 
% * input lags $x_u$ with length $l_u$
% * vector $y$ of single output,
% * output lags $x_y$  with length $l_y$
% * with FCM membership function
% $$ \mu_i(z) = \left( \sum_{j=1}^{n_v} \left( \frac{||z-v_i||}{ ||z-v_j||} \right)^{\dfrac{2}{\nu-1}} \right)^{-1} $$
% * or Gauss membership function
% $$ \mu_i(z) = e^{-\dfrac{||z-v_i||^2}{2\cdot\sigma_i^2}} $$
% * norm $$ ||z-v_j|| = (z-v_j)^T\cdot w_j\cdot (z-v_j) $$
% * and fuzzy basis functions
% $$ \phi_i(z) = \frac{\mu_i(z)}{\sum_{j=1}^{n_v} \mu_j(z)} $$
% * with the scheduling variable $z=u$ (for input space clustering) or $z=[u,y]$ (for product
% space  clustering), and
% * cluster centers $v_i, i=1,\ldots,n_v$.

%% Algorithm
%%
% # Choose number of local models: $n_v=3$  and fuzzy parameter $\nu=1.2$ or $\sigma_i=???$
% # Select the TS model with minimal MSE of $s$ multi-start tries with clustering and
% global LS-estimation. 
% # Optimize the TS model parameters $(v_i,B_i,c_i)$ for each try or the best found model.

%% Structural parmeters
nv = 3;    % number of local models    
nu = 1;    % number of inputs   
ny = 1;    % number of outputs  
nue = 1.2; % Fuzziness parameter
%%
% Scheduling lags $z_{lag}$ = regressor lags $x_{lag}=[y_{k-1},y_{k-2},u_k]^\top$ 
z_lag_u = {0};
z_lag_y = [1,2];
x_lag_u = {0};
x_lag_y = [1,2];

%% Identification data
%
% Create input $u$ as steps with amplitude 2 and random 
% width $l=[1,\ldots,20]$ for $N=1000$ time
% steps (sampling rate is 0.01s) and compute the output $y$ from the
% Narendra function 
N = 1000;
dt = 1e-2;                    % Sampling time
t = dt * transpose( 0:N-1 );  % time vector: $t$
rng(0);                       % Initalize random number generator
[ u, y ] = Narendra_fct( N );
%%
% Plot of the identification data
figure(1);clf

subplot(2,1,1)
plot(t,u)
grid on
ylabel('u')
title(sprintf('Narendra function / input: steps / N=%d',N))

subplot(2,1,2)
plot(t,y)
grid on
ylabel('y')
xlabel('t')

set(gcf,'WindowState', 'maximized' );

%% Create the NARX TS model
%
addpath( '../TSModel' );  % Path to TSModel class
ts = TSModel( 'ARX', nv, nu, 'Name','NARX Narendra', ...
    'Comment','Narendra function');
%%
% Set the identification data $u$, $y$
ts.setData( u, y, 'SampleTime',dt, 'Labels', { 'u', 'y' } );
%%
% Set the data limits: u=[-2,2], y=[-5,10]
ts.setDataLimits( [-2,2 ; -5,10] );
%%
% Set the scheduling and regressor lags 
ts.setSchedulingLags( z_lag_u, z_lag_y );
ts.setRegressorLags( x_lag_u, x_lag_y );

%% Clustering 
% Clustering in product-space $z=[u,y]$ with FCM membership functions
ts.clustering( 'FCM', 'nue',nue, 'tries',1, 'seed',0 );
%%
% Estimated cluster centers: colums $y_{k-1},y_{k-1},u_k$
v1 = getCluster( ts )

%% Initialize the local models
% Estimate the local model parameters with global Least Squares
ts.initialize( 'FCM', 'nue', nue, 'method','global'  );

%% Predict the  TS model output
% for given data $u$: $y_{pred}$
y_pred = ts.predict( u, y );
%%
% Plot the correlation
plotResiduals( y, y_pred, 'figure', 2, ...
    'title', 'Narendra/NARX: correlation ' );
set(gcf,'WindowState', 'maximized' );
%% 
% Plot of identification data $y$ and predicted data $y_{pred}$
figure(3);clf

yyaxis left
plot(t,u,'b--')
ylabel('u')

yyaxis right
plot(t,y,'g-',t,y_pred,'r-')

grid on
ylabel('y')
xlabel('t')
title( sprintf('Narendra/NARX: n_v=%d, \\nu=%g',ts.nv,ts.nue))
legend('u','y_{obsv}','y_{pred}')
set(gcf,'WindowState', 'maximized' );

%% Validate with new data
% Test with new step inputs $u_{val}$, $y_{val}$
[u_val,y_val] = Narendra_fct( N );
y_pred_val = ts.predict( u_val,y_val );
%%
% Plot of observed and predicted outputs
figure(5),clf

yyaxis left
plot(t,u_val,'b--')
ylabel('u')

yyaxis right
plot(t,y_val,'g-',t,y_pred_val,'r-')

grid on
ylabel('y')
xlabel('t')
title( sprintf('Narendra/NARX validation: n_v=%d, \\nu=%g',ts.nv,ts.nue))
legend('u','y_{obsv}','y_{pred}')
set(gcf,'WindowState', 'maximized' );
%%
% Plot the correlation
plotResiduals( y, y_pred_val, 'figure', 4, 'title', 'Narendra/NARX: correlation' );
set(gcf,'WindowState', 'maximized' );

%% Optimize the TS model parameters 
% Optimize both, clusters centers $v$ (MF) and local model parameters $A_i,B_i,c_i$
ts.optimize( 'B' )
%%
% Cluster centers of the optimized TS model
v2 = getCluster( ts )
%%
% Plot the correlation
y_pred_opt = ts.predict( u,y );
plotResiduals( y, y_pred_opt, 'figure', 5, ...
    'title', 'Narendra/NARX opt: correlation' );
set(gcf,'WindowState', 'maximized' );
%% 
% Plot of observed and predicted outputs for optimized TS model
figure(6);clf

yyaxis left
plot(t,u,'b--')
grid on
ylabel('u')
title( sprintf('Narendra/NHARX: optimized model n_v=%d \\nu=%g',ts.nv,ts.nue))

yyaxis right
plot(t,y,'g-',t,y_pred_opt,'r-')
grid on
ylabel('y')
xlabel('t')
legend('u','y_{obsv}','y_{pred}')
set(gcf,'WindowState', 'maximized' );
