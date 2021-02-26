%% Takagi-Sugeno Model Identification Toolbox
%
% Example of a NARX LiP TS model for the Narendra function
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
% Example of the identification of a NARX MISO LiP TS model for 
% given inputs $u$ and output $y$.
%% 
% Determine the NARX LiP TS model
%
% $$ y(k+1) = \sum_{i=1}^{n_v} \phi_i(z) \cdot \left( \sum_{l=0}^{l_y} A_{i,j}\cdot y_{k-l} + \sum_{j=0}^{n_u}\sum_{j=0}^{l_u} B_{i,j}\cdot u_{k-l} + c_{i} \right)$$
%%
% * for given $u_j, j=1,\ldots,n_u$ of $n_u$ input vectors and 
% * input lags $l_u$
% * vector $y$ of single output,
% * output lags $l_y$
% * with FCM membership function
% $$ \mu_i(z) = \left( \sum_{j=1}^{n_v} \left( \frac{||z-v_i||}{ ||z-v_j||} \right)^{\dfrac{2}{\nu-1}} \right)^{-1} $$
% * or Gauss membership function
% $$ \mu_i(z) = e^{-\dfrac{||z-v_i||^2}{2\cdot\sigma_i^2}} $$
% * norm ||z-v_j|| = (z-v_j)^T\cdot w_j\cdot (z-v_j)
% * and fuzzy basis functions
% $$ \phi_i(z) = \frac{\mu_i(z)}{\sum_{j=1}^{n_v} \mu_j(z)} $$
% * with the scheduling variable $z=u$ (for input space clustering) or $z=[u,y]$ (for product
% space  clustering), and
% * cluster centers $v_i, i=1,\ldots,n_v$.
%%
% Algorithm:
%%
% # Select the TS model with minimal MSE of $m$ multi-start tries with clustering and
% LS-estimation. 
% # If no $n_v$ is choosen: try for $n_v=3,\ldots,n{v,\max}$
% # If no fuzzy parameter is choosen, try $\nu=[1.05,1.1,1.2,1.5,2]$ or $\sigma_i=???$
% # Optimize the TS model parameters $(v_i,B_i,c_i)$ for each try or the best found model.

%% Approximation of the Narendra function
% 
% $$ y_{k+1} = \frac{ y_k\cdot y_{k-1} \cdot ( y_k + 2.5 ) }{ 1 + y_k^2 + y_{k-1}^2 } +  u_k$$
%
% with inital states $y_1 = 0, y_2 = 0$
%
% as a NARX TS model with $n_v=3$ local models
%%
% Source; K. Narendra and K. Parthasarathy. "Identification and control of dynamical systems using neural networks". IEEE Transactions on Neural Networks 1(1) (Mar. 1990), pp. 4-27.
%
% <https://www.sciencedirect.com/science/article/pii/0888613X9290014Q>


%% Structural parmeters
nc = 3;    % number of local models    
nu = 1;    % number of inputs   
ny = 1;    % number of outputs  
nue = 1.2; % Fuzziness parameter
%%
% Scheduling lags $z_{lag}$ = regressor lags $x_{lag}=[y(k-1),y(k-2),u(k)]$ 
z_lag_u = {0};
z_lag_y = [1,2];
x_lag_u = {0};
x_lag_y = [1,2];

%% Compute identification data
%
% Create input $u$ as steps with width $l=1,\ldots,20]$ for $N=1000$ time steps 
% and compute the output $y$ from the Narendra function 
N = 1000;
rng(0); % Initalize random number generator
[ u, y ] = Narendra_fct( N );

%%
% Sampling time
dt = 1e-2; 
%%
% Compute the time vector: $t$
t = dt * transpose( 0:size(u,1)-1 ); 

%% Plot of the identification data
h=figure(1);clf

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

h.WindowState = 'maximized';

%% Creation of the NARX TS model
ts = TSModel( 'ARX', nc, nu, 'Name','NARX Narendra', 'Comment','Narendra function');
%%
% Set the scheduling and regressor lags 
ts.setSchedulingLags( z_lag_u, z_lag_y );
ts.setRegressorLags( x_lag_u, x_lag_y );
%%
% Set the identification data $u$,$y$
ts.setData( u, y, 'SampleTime',dt, 'Labels', { 'u', 'y' } );
%%
% Set the data limits: u=[-2,2], y=[-5,10]
ts.setDataLimits( [-2,2 ; -5,10] );

%% Clustering in product-space: $z=[u,y]$ with FCM membership functions:
ts.clustering( 'FCM', 'nue',nue, 'tries',1, 'seed',0 )

%% Get the found cluster centers
v1 = getCluster(ts)

%% Initialize the local models with global LS:
ts.initialize( 'FCM', 'nue', nue, 'method','global'  );

%% Predict the  TS model ouput for given data: $y_{pred}$
y_pred = ts.predict( u, y );
%%
% Plot the residuals
hr=plotResiduals( y, y_pred, 'figure', 2, 'title', 'Correlation Narendra/NARX' );
hr.WindowState = 'maximized';

%% 
% Plot identification data $y$ and predicted data $y_{pred}$
h=figure(3);clf

yyaxis left
plot(t,u,'b--')
ylabel('u')

yyaxis right
plot(t,y,'g-',t,y_pred,'r-')

grid on
ylabel('y')
xlabel('t')
title( sprintf('Narendra: NARX model n_v=%d \\nu=%g',ts.nv,ts.nue))
legend('u','y_{obsv}','y_{pred}')
h.WindowState = 'maximized';

%% Predict output from new test data $u_t$, $y_t$
[ut,yt] = Narendra_fct( N );
yp_t = ts.predict( ut,yt );

%%
% Plot the correlation
hc=plotResiduals( y, yp_t, 'figure', 4, 'title', 'Residuals Narendra eval/NARX' );
h.WindowState = 'maximized';
%%
% Plot 
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
h.WindowState = 'maximized';

%% Optimization of TS model parameters 
% Optimize clusters centers $v$ (MF) and local model parameters $A_i,B_i,c_i$
ts.optimize( 'B' )
%%
% Get cluster centers of optimized TS model
v2 = getCluster( ts )
%%
% Plot the residuals
yp_o = ts.predict( u,y );
hr = plotResiduals( y, yp_o, 'figure', 6, 'title', 'Residuals Narendra/NARX opt' );
hr.WindowState = 'maximized';

%% 
% Plot identification data and predicted output of optimized TS model
h=figure(7);clf

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
h.WindowState = 'maximized';
