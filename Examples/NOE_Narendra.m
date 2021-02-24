%% Takagi-Sugeno Model Identification Toolbox
%
% NOE LiP model example for the Narendra function
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

% $Id$%


addpath( '..' ) % Path to TSModel files

%% Model order
nc = 3;    % number of clusters = local models
nu = 1;    % number of inputs
ny = 1;    % number of outputs
nue = 1.2; % fuzziness parameter

%% Generate Narendra model data
n = 1000;
rng(0);
[u,y] = Narendra_fct( n );

dt = 1e-2; % Sampling time
t = dt * transpose( 0:size(u,1)-1 );

figure(1),clf
subplot(2,1,1)
plot(t,u)
grid
ylabel('u')
subplot(2,1,2)
plot(t,y)
grid
ylabel('y')
xlabel('t')

%% Create TS model
ts = TSModel( 'OE', nc, nu, 'Name','OE Narendra', 'Comment','Narendra function');
ts.setLags( [0], [2] );

ts.setData( u, y, 'SampleTime',dt, 'Labels', { 'u', 'y' } );
ts.setDataLimits( [-2,2 ; -5,10] );

%% Clustering in product-space (u,y)
ts.clustering( 'FCM', 'nue', nue, 'tries',1, 'seed', 0 )
v1 = getCluster(ts);

%% Initialize local LS models with LS (global/local)
ts.initialize( 'FCM', 'nue', nue, 'method','global'  );
%ts.initialize( 'FBF', 'nue', nue, 'method','local'  );

%% Compute TS model for given data
yp = ts.predict( u,y );
plotResiduals( y, yp, 'figure', 2, 'title', 'Narendra NARX: correlation' );

%%
figure(3),clf
plot(t,u,'k-',t,y,'g-',t,yp,'r-')
grid on
title('Narendra NOE: observed data')
legend('u','y_{obsv}','y_{pred}')

%% Eval on unkown data
[ut,yt] = Narendra_fct( n );

ypt = ts.predict( ut,yt );
plotResiduals( y, ypt, 'figure', 4, 'title', 'Narendra NOE: corrleation on test data' );

figure(5),clf

yyaxis left
plot(t,ut,'b--')
ylabel( 'u' ) 
yyaxis right

plot(t,yt,'g-',t,ypt,'r-')
ylabel( 'y' ) 
xlabel( 't' )

grid on
title('Narendra NOE: eval data')
legend('u','y','y_{pred}')

%% Optimize Clusters c (MF) and/or local model A/B/C
ts.optimize( 'B' )
v2 = getCluster( ts );

ypo = ts.predict( u,y );
plotResiduals( y, ypo, 'figure', 6, 'title', 'Narendra NOE: correlation opt' );

figure(7),clf
plot(t,ut,'k-',t,yt,'g-',t,ypt,'r-')
grid on
title('Narendra NOE: opt predicted ouput')
legend('u','y','y_{pred}')

