%% Takagi-Sugeno Model Identification Toolbox
%
% Automatic static LiP model for an academic example
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

%%
% 
% Example of automatic identification of a static MISO LiP TS model for 
% given multiple inputs $u$ and single output $y$ with minimal requirements.

%% 
% Determine the MISO LiP TS model
%
% $$y(u) = \sum_{i=1}^{n_v} \phi_i(z) \cdot \left(\sum_{j=1}^{n_u} B_{i,j}\cdot u_j + c_{i}\right) $$
%%
% * for given vectors $u_j, j=1,\ldots,n_u$ of $n_u$ inputs and 
% * vector $y$ of single output,
% * with FCM membership function
% $$ \mu_i(x) = \left( \sum_{j=1}^{n_v} \left( \frac{||z-v_i||}{ ||z-v_j||} \right)^{\dfrac{2}{\nu-1}} \right)^{-1} $$
% * or Gaussian membership function
% $$ \mu_i(z) = e^{-\dfrac{||z-v_i||^2}{2\cdot\sigma_i^2}} $$
% * norm $||z-v_j|| = (z-v_j)^T\cdot A_j\cdot (z-v_j)$
% * and fuzzy basis functions
% $$ \phi_i(z) = \frac{\mu_i(z)}{\sum_{j=1}^{n_v} \mu_j(z)} $$
% * with the scheduling variable $z=u$ (for input space clustering) or $z=[u,y]$ (for product
% space  clustering), and
% * cluster centers $v_i, i=1,\ldots,n_v$.

%% Algorithm
%%
% # Search the best TS model with the minimal MSE for $n_v=\{2,3,4\}$ and $\nu=\{1.05,1.1,1.2,1.5,2\}$.
% # Select the TS model with minimal MSE of $s$ multi-start tries for clustering and
% Least Squares estimation. 
% # Optimize the TS model parameters $(v_i,B_i,c_i)$ for each try.

%% Minimal required data
%
% Inputs $u\in\mathbb{R}^{N \times n_u}$ and output $y\in\mathbb{R}^N$, each with $N$ data points

%% Identification data 
%
% Given is an academic example as a TS model with 
%% 
% * inputs $u_1,u_2\in[0,2]$ 
% * local model matrices 
% $$ B = \begin{pmatrix} -4  & 4 \\ 4 & -2 \\ 2 & 1\end{pmatrix},\quad c = \begin{pmatrix} -2\\-4\\1 \end{pmatrix} $$ 
% * FCM membership functions ($\nu=1.2$) with Euclidean norm
% * cluster centers 
% $$ v = \begin{pmatrix} 0.5 & 0.5 \\ 0.5 & 1.5 \\ 1.5 & 1 \end{pmatrix} $$

%%
% Load data $u$, $y$ with $N=50$ data-points without noise, generated from this model: 
load( 'Data/AcadEx.mat' )

%% Structural parameters
% Number of inputs $n_u$ = number of columns in $u$
Par.nu = size( u, 2);    
%%
% Number of clusters $n_v$ = number of local models ($n_v$ > 1)
Par.nv = [ 2, 3, 4 ];
%%
% Fuzziness parameter (FCM: $\nu = \{1.05,\ldots, 2\}$, Gauss: $\sigma_i^2$)
Par.fuzzy = [ 1.05, 1.2 ,2.0 ];

%% Optional settings
% 
% For more control over the approximation process.
%% 
% Multi-Start: number of tries $s$ (clustering & LS), default = 10
Par.Tries = 10;
%%
% Clustering: Fuzzy C-Means (FCM) / Gustafson-Kessel (GK) / KMeans (KMeans), default = 'FCM'
Par.Clustering = 'FCM';
%%
% Clustering in product space:  $u$ and $y$ (true) or only input space $u$ (false)
Par.ProductSpace = true;
%%
% Norm for clustering: 'Euclidean' or 'Mahalanobis', default = 'Euclidean'
Par.Norm = 'Euclidean';
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
% Debug infos of algorithm progess: (0=none, 1=info, 2=detailed)
Par.Debug = 2;

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

%% Retrieve the parameters of the final TS model 
%%
% Show the TS model parameters:
disp( model )
%%
% Show the cluster centers $v$ ($n_v$ rows and $n_u$ columns): 
v = getCluster( model )
%%
% Show the local model matrices $B_i$ and $c_i$ ($n_v$ rows and $n_u$ columns):
[~,B,c] = getLM( model )