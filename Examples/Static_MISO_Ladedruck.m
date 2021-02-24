%% TS-Toolbox: LiP MISO example (Ladedruck )

% $Id$

clear
close all

addpath( '../TSModel' ) % Path to TSModel files

%% Modell order
nc = 3;    % number of clusters = local models
% Selection of inputs in data
inu = [1,2,3,4,5]
%inu = 1:5;
nu = length(inu);    % number of inputs
nue = 1.2; % fuziness parameter

% Input vector u and output vector y
load( 'Data/MISO_Ladedruck.mat' )
dt = 1; % Implicit sampling time for static models

%% Create TS model
ts = TSModel( 'ARX', nc, nu, 'Name','Ladedruck',...
'Comment','IAV Überlandfahrt 24.3.2013 Konstante Einspritzcharakteristik');
ts.setFuziness( nue );
ts.set_msf_type( 'FBF' );
ts.setSchedulingLags( {0,0,0,0,0}, [1,2,3] );
ts.setRegressorLags( {0,0,0,0,0}, [1,2,3] );

%% I/O data
Labels = {};
u = [];
for i=inu
    Labels{end+1} = data.u(i).label;
    u = [u,data.u(i).Vals];
end
Labels{end+1} = data.y.label;
y = data.y.Vals;
dt = data.y.ts;
ts.setData( u, y, 'SampleTime',dt, 'Labels', Labels,  ...
    'Comment', 'from dataset' );
%ts.setDataLimits( [0,2 ;0,2; -7,7] );

%% Clustering in product-space (u,y)
ts.clustering( 'FCM', 'nue', nue, 'tries',2, 'seed', 0 )
c1 = getCluster(ts);

%% Plot co Clustering
figure(1),clf
x = [y,u];
c = c1(:,[1,4:end]);
[sr,sc] = getSubplotPar( nu + 1 );
s=0;
for i1=1:nu
    for i2=i1+1:nu+1
        s=s+1;
        subplot(sr,sc,s)
        plot(x(:,i1),x(:,i2),'.')
        grid on
        xlabel( ts.Labels(i1) )
        ylabel( ts.Labels(i2) )
        hold on
        plot(c(:,i1),c(:,i2),'x')
        hold off
    end
end
orient landscape
print('-dpdf','-bestfit','Test_SISO_Ladedruck_clustering.pdf')

%% Initial modell
ts.initialize( 'FBF', 'nue', nue, 'method','global'  );

%% Evaluation of initial modell
yp = ts.evaluate( u,y );
plotResiduals( y, yp, 'figure', 2, 'title', 'Residuals Throttle NARX' );
orient landscape
print('-dpdf','-bestfit','Test_SISO_Ladedruck.pdf')

return

%% Optimize Clusters c (MF) and/or local model A/B/C
ts.optimize( 'B' )
c2 = getCluster( ts );
ypo = ts.evaluate( u,y );
plotResiduals( y, ypo, 'figure', 3, 'title', 'Residuals MISO Ladedruck NARX opt' );

orient landscape
print('-dpdf','-bestfit','Test_SISO_Ladedruck_opt.pdf')