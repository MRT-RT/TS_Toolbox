%% Takagi-Sugeno Model Identification Toolbox
%
% Static LiP model for the 2-dimensional-Friedman function.
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

%% Identification data 
% Use the 2D-Friedman function:
% $$ y = 10\cdot\sin( \pi\cdot u_1 \cdot u_2 ) $$
nu = 2;
%%
% Choose input matrix $u$ as random data with $N$ data-points: $u_{1,2}\in[0,1]$
N = 500;
u = rand( N, nu );
%%
% Compute output vector $y$ from the Friedman function:
y = Friedman_fct( u, nu );

%% Structural parameters
% Number of inputs $n_u$ = number of columns in $u$
Par.nu = size( u, 2);    
%%
% Number of clusters $n_v$ = number of local models ($n_v$ > 1)
Par.nv = 0;    
%%
% Fuzziness parameter (FBF: $\nu = [1.05,\ldots, 2]$, Gauss: $\sigma^2$)
Par.fuzzy = 1.2; 

%% Optional settings
% 
% For more control over the approximation process.
%% 
% Multi-Start: number of tries $m$ (clustering & LS), default = 10
Par.Tries = 3;
%%
% Clustering: Fuzzy C-Means (FCM) / Gustafson-Kessel (GK) / KMeans (KMeans), default = 'FCM'
Par.Clustering = 'FCM';
%%
% Clustering in product space:  $u$ and $y$ (true) or only input space $u$ (false)
Par.ProductSpace = true;
%%
% Norm for clustering: 'Euclidian' or 'Mahalanobis', default = 'Euclidian'
Par.Norm = 'Euclidan';
%%
% Membership functions: 'FCM' clustering or 'Gauss' type
Par.MSF = 'FCM';
%%
% Least Squares estimation of local models: 'local' or 'global', default = 'global'
Par.LS = 'global';
%%
% Optimize model parameters: default='both'
%%
% * no optimization: 'none',
% * only $v$: 'cluster',
% * only local models ($B_i,c_i$): 'model', or
% * both $v$ and $B_i,c_i$: 'both'
Par.Optimize = 'both';
%%
% Optimize each try or only best try: default='each'
%%
% * each try: 'each',
% * best try: 'best' (less computation time)
Par.IterOpt = 'each';
%%
% Plot clusters and residuals: 'none'/'iter'/'final', default='final'
Par.Plots = 'final';

% Debug infos (0=none, 1=info, 2=detailed)
Par.Debug = 2;

%% Estimation of  LiP TS model parameters
%%
% Estimate the TS model with plot of clustering and correlation:
model = TSM_Static_auto( u, y, Par );
%%
% Predict the model output $y_{pred}$ for input $u$:
y_pred = model.predict( u, y );
%%
% Plot a residual histogram:
hr = plotResidualHist( y, y_pred, 'figure', 3 );
%%
% Plot the rule activation and input/output data:
ha = plotRuleActivation( u,y,model, 'figure', 4 );
%%
% Show the parameter of the resulting TS model:
disp( model )

%% Test with unknown data
u_test = rand( N, nu );
y_obsv = Friedman_fct( u_test, nu );%%
% Calculate output vector $y_{pred}$
y_pred = model.predict( u_test );

figure(3),clf
plot( 1:N, y_obsv, 'k-',1:N, y_pred, 'r--' )
grid on
xlabel('k')
ylabel('y')
title( 'Friedman-2D: observed vs. predicted ouput' )
legend( 'y_{obsv}','y_{pred}' )

%%
h = plotResiduals( y_obsv, y_pred, 'figure', 4, 'title', 'Friedman-2D: correlation' );
