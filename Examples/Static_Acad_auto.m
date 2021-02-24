%% Takagi-Sugeno Model Identification Toolbox 
%
% Automatic static LiP model for an academicx example
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
% 
% Example of automatic identification of a static MISO LiP TS model for 
% given inputs $u$ and output $y$ with minimal requirements.

%% 
% Determine the MISO LiP TS model
%
% $$y(u) = \sum_{i=1}^{n_v} \phi_i(z) \cdot \left(\sum_{j=1}^{n_u} B_{i,j}\cdot u_j\right) + c_{i}$$
%%
% * for given $u_j, j=1,\ldots,n_u$ of $n_u$ input vectors and 
% * vector $y$ of single output,
% * with FCM membership function
% $$ \mu_i(x) = \left( \sum_{j=1}^{n_v} \left( \frac{||z-v_i||}{ ||z-v_j||} \right)^{\dfrac{2}{\nu-1}} \right)^{-1} $$
% * or Gauss membership function
% $$ \mu_i(z) = e^{-\dfrac{||z-v_i||^2}{2\cdot\sigma_i^2}} $$
% * norm ||z-v_j|| = (z-v_j)^T\cdot A_j\cdot (z-v_j)
% * and fuzzy basis function
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

%% Minimal required data
%
% * Inputs $u\in\mathbb{R}^{N \times n_u}$ and output $y\in\mathbb{R}^N$, each with $N$ data points

%% Additional choices
%
% # Number of local models $n_v > 2$, recommended initial choice $n_v=3$
% # Fuzziness parameter $\nu = [1.05,\ldots,2]$, recommended initial choice: $\nu=1.2$

%% Identification data 
% Given is the academic example as a TS model with 
%% 
% * inputs $u_1,u_2\in[0,2]$ 
% * local model matrices $B,c$ 
% $$ B = \begin{pmatrix} -4  & 4 \\ 4 & -2 \\ 2 & 1\end{pmatrix},\quad c = \begin{pmatrix} -2\\-4\\1 \end{pmatrix} $$ 
% * FCM membership functions ($\nu=1.2$)
% * cluster centers 
% $$ v = \begin{pmatrix} 0.5 & 0.5 \\ 0.5 & 1.5 \\ 1.5 & 1 \end{pmatrix} $$

%%
% Load data $u$, $y$ generated with this model from file: 
load( 'Data/AcadEx.mat' )

%% Structural parameters
% Number of inputs $n_u$ = number of columns in $u$
Par.nu = size( u, 2);    
%%
% Number of clusters $n_v$ = number of local models ($n_v$ > 1)
Par.nv = 3;     % choose a fixed value
% Par.nv = 2:5; % choose a range "from : to"
% Par.nv = 0;   % autoselect range

%%
% Fuzziness parameter (FCM: $\nu = [1.05,\ldots, 2]$, Gauss: $\sigma_i^2$, 0=select range)
Par.fuzzy = 1.2;                   % choose a fixed value
%Par.fuzzy = [1.05, 1.2, 1.6, 2];  % choose set of values
%Par.fuzzy = 0;                    % autoselect range

%% Optional settings
% 
% For more control over the approximation process.
%% 
% Multi-Start: number of tries $m$ (clustering & LS), default = 10
Par.Tries = 10;
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
% Membership functions: 'FCM' or 'Gauss' type clustering
Par.MSF = 'FCM';
%%
% Least Squares estimation of local models: 'local' or 'global', default = 'global'
Par.LS = 'global';
%%
% Optimize TS model parameters: default='both'
%%
% * no optimization: 'none',
% * only $v$: 'cluster',
% * only local models ($B_i,c_i$): 'model', or
% * both $v$ and $B_i,c_i$: 'both'
Par.ParOpt = 'both';
%%
% Optimize each try or only best try: default='each'
%%
% * each try: 'each',
% * best try: 'best' (less computation time)
Par.IterOpt = 'each';
%%
% Plot clusters and residuals: 'none'/'iter'/'final', default='final'
Par.Plots = 'final';
%%
% Debug infos (0=none, 1=info, 2=detailed)
Par.Debug = 1;

%% Estimation of  Static TS model parameters
%%
% Estimate the TS model with plot of clustering and correlation:
model = TSM_Static_auto( u, y, Par );
%%
% Predict the model output $y_{\mathrm{pred}}$ for input $u$:
y_pred = model.predict( u, y );
%%
% Plot a residual histogram:
hr = plotResidualHist( y, y_pred, 'figure', 3 );
%%
% Plot the rule activation and input/output data:
ha = plotRuleActivation( u,y,model, 'figure', 4 );

%% Retrieve the TS model parameter
%%
% Show the TS model parameters:
disp( model )
%%
% Show the cluster centers $v$ ($n_v$ rows and $n_u$ columns): 
v = getCluster( model )
%%
% Show the local model matrices $B_i$ and $c_i$ ($n_v$ rows and $n_u$ columns):
[~,B,c] = getLM( model )

%% Prediction of the TS model 
%
%%
% As test data, use random inputs $u_{test}\in [0,2]\times[0,2]$ 
% with $N=2000$ data points 
u_test = 2 * rand( 2000, Par.nu );
y_test = model.predict( u_test );
%%
% Plot of model outputs for random test input data:
he = figure( 10 );
plot3( u(:,1), u(:,2),y,'b.', ...
       u_test(:,1), u_test(:,2), y_test,'r.','MarkerSize',8);
legend( 'u_{ident}', 'u_{test}','Location','NW' )
grid on
xlabel('u_1'),ylabel('u_2'),zlabel('y')
title( 'Static auto: random test data' )
set(gca,'FontSize', 14)
he.WindowState = 'maximized';

%%
% For grid type distributed test data: use $n$ points per input $u_i$
n = 40;
[U1,U2] = meshgrid( linspace(0,2,n), linspace(0,2,n) );
u_grid = [ U1(1:end)', U2(1:end)' ];
y_grid = model.predict( u_grid );
Y = reshape( y_grid,n,n );
%%
% Plot of model output for grid input data:
hg = figure( 11 );
plot3( U1, U2, Y, 'k.')
grid on
xlabel('u_1'),ylabel('u_2'),zlabel('y')
title( 'Static auto: grid data' )
set(gca,'FontSize', 14)
hg.WindowState = 'maximized';

%% Plot the final TS model in 2D
% Get the membership degrees for the grid data
mu = getMSF( model, u_grid, y_grid );

figure(12),clf
hold on
plot( v(:,1), v(:,2),'rx' )
for i=1:model.nv
    contour( U1,U2, reshape(mu(:,i),n,n), 1 )
end
hold off
axis square
grid on, box on
view(0,90)
xlabel('u_1'),ylabel('u_2'),zlabel('\mu(z)')
title( 'Cluster centers' )
legend( 'v', '\mu_i=0.5' )

set(gca,'FontSize', 14)
hg.WindowState = 'maximized';

figure(13),clf
hold on
for i=1:model.nv
    meshc( U1,U2, reshape(mu(:,i),n,n) )
end
hold off
axis square
view(-20,40)
xlabel('u_1'),ylabel('u_2'),zlabel('\mu(z)')
title( 'Membership degrees' )
set(gca,'FontSize', 14)
set(gca,'FontSize', 14)
hg.WindowState = 'maximized';

figure(14),clf
hold on
plot3( u(:,1),u(:,2),y,'k.' )
surf( U1,U2, reshape(Y,n,n) )
hold off
title( 'Predicted TS model' )
xlabel('u_1'),ylabel('u_2'),zlabel('y')
grid on, box on
axis square
view(-20,40)
xlabel('u_1'),ylabel('u_2'),zlabel('\mu(z)')
%legend( 'ident data','prediction' )

set(gca,'FontSize', 14)
hg.WindowState = 'maximized';

