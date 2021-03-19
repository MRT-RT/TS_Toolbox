%% Takagi-Sugeno Model Identification Toolbox
%
% Example of a NOE TS model for the Narendra function
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
% $Id$%

%%
% Example of the identification of a NOE MISO TS model for 
% given multiple inputs $u$ and single output $y$.
%% 
% Determine the NOE TS model
%
% $$ \hat{y}_{k+1} = \sum_{i=1}^{n_v} \phi_i(z) \cdot \left( \sum_{l=0}^{l_y} A_{i}\cdot \hat{y}_{k-l} + \sum_{j=0}^{n_u}\sum_{l=0}^{l_u} B_{i,j}\cdot u_{k-l} + c_{i} \right) + e_k$$
%%
% * for given $u_j, j=1,\ldots,n_u$ of $n_u$ input vectors, output error $e$ and 
% * input lags $x_u$ with length $l_u$
% * intial vector $y$ of single output,
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


%% Structural settings
nc = 3;    % number of clusters = local models
nu = 1;    % number of inputs
ny = 1;    % number of outputs
nue = 1.2; % fuzziness parameter
%%
% Scheduling lags $z$ = regressor lags $x=[\hat{y}_{k-1},\hat{y}_{k-2},u_k]^\top$ 
z_lag_u = {0};
z_lag_y = [2];
x_lag_u = {0};
x_lag_y = [2];

%% Identification data
%
% Create input $u$ as steps with width $l=[1,\ldots,20]$ for $N=1000$ time
% steps (sampling rate is 0.01) and compute the output $y$ from the Narendra function 
N = 1000;
dt = 1e-2;                           % Sampling time
t = dt * transpose( 0:size(u,1)-1 ); % time vector $t$
rng(0);
[u,y] = Narendra_fct( N );
%% 
% Plot of the identification data
h=figure(1);clf

subplot(2,1,1)
plot(t,u)
grid on
ylabel('u')
subplot(2,1,2)
plot(t,y)
grid  on
ylabel('y')
xlabel('t')

%% Creation of  TS model
addpath( '../TSModel' );  % Path to TSModel class
ts = TSModel( 'OE', nc, nu, 'Name','OE Narendra', 'Comment','Narendra function');
ts.setSchedulingLags( z_lag_u, z_lag_y );
ts.setRegressorLags( x_lag_u, x_lag_y );
%%
% Set the identification data
ts.setData( u, y, 'SampleTime',dt, 'Labels', { 'u', 'y' } );
ts.setDataLimits( [-2,2 ; -5,10] );

%% Clustering 
% Clustering in product-space $z=[u,y]$ with FCM membership functions and
% $\nu=1.2$ with $s=3$ multi-start tries and ficed initialized random number
% generator (seed 0)
ts.clustering( 'FCM', 'nue', nue, 'tries',3, 'seed', 0 )
%%
% Cluster centers of the inital model
v1 = getCluster( ts )

%% Initialization of local models 
% with global Least-Squares, FCM membership functions and $\nu=1.2$
ts.initialize( 'FCM', 'nue', nue, 'method','global'  );

%% Predicted NOE TS model ouput
y_pred = ts.predict( u,y );
plotResiduals( y, y_pred, 'figure', 2, 'title', 'Narendra NOE: correlation' );
set(gcf,'WindowState', 'maximized' );

%%
% Plot of the observed vs. predicted outputs
figure(3);clf

plot(t,u,'k-',t,y,'g-',t,y_pred,'r-')
grid on
xlabel( 'time t' )
ylabel( 'output y' )
title('Narendra NOE: observed vs. predicted outputs')
legend('u','y_{obsv}','y_{pred}')
set(gcf,'WindowState', 'maximized' );

%% Prediction on validation data
[u_val,y_val] = Narendra_fct( N );
y_val_pred = ts.predict( u_val,y_val );
%%
% Plot the correlation
plotResiduals( y, y_val_pred, 'figure', 4, 'title', 'Narendra NOE: corrleation on test data' );
set(gcf,'WindowState', 'maximized' );
%%
% Plot of observed vs. predicted validation data
figure(5);clf

yyaxis left
plot(t,u_val,'b--')
ylabel( 'u' ) 
yyaxis right

plot(t,y_val,'g-',t,y_val_pred,'r-')
ylabel( 'y' ) 
xlabel( 't' )

grid on
title('Narendra NOE: validation data')
legend('u','y_{obsv}','y_{pred}')
set(gcf,'WindowState', 'maximized' );

%% Optimize the TS model parameters
%
% Set additional parametrs for function |lsqnonlin| 
optimopts = optimoptions('lsqnonlin');
optimopts.FunctionTolerance = 1e-6;
optim.OptimalityTolerance = 1e-6;
optimopts.StepTolerance = 1e-12;
optimopts.Display = 'iter-detailed';
%%
% Optimize both, the cluster centers $v$ (MF) and 
% the local model parameters $A_i, B_i, c_i$
ts.optimize( 'Both', 'optimopts', optimopts );
%%
% Get the cluster centers of the optimized NOE TS model
v2 = getCluster( ts )
%%
% Plot the correlation on the validation data
y_pred_opt = ts.predict( u,y );
plotResiduals( y, y_pred_opt, 'figure', 5, 'title', 'Narendra NOE: correlation opt' );
set(gcf,'WindowState', 'maximized' );
%%
% Plot the observed vs. the predicted validation data
figure(6);clf
yyaxis left
plot(t,u_val)
ylabel( 'input u' )
yyaxis right
plot(t,y_val,'g-',t,y_val_pred,'r-')
grid on
ylabel( 'output y' )
xlabel( 'time t' )
title('Narendra NOE opt: predicted validation output')
legend('u','y_{obsv}','y_{pred}')
set(gcf,'WindowState', 'maximized' );

