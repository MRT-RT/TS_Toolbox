%% plotMSF_activation
%
% $Id$

function plotMSF_activation( tsm, u )

[n,nu] = size(u);
    
mu = tsm.z_msf( u, tsm.c, tsm.m );

figure(1)
clf
xlim([0,n+1])
ylim([-0.05,1.05])
col = colormap( jet(tsm.nc) );
l = zeros(n,tsm.nc);
for i = 1 : n
    l(i,1) = mu(i,1);
    line( [i,i], [0,l(i,1)],'color', col(1,:) )
    for j=2 : tsm.nc
        l(i,j) = l(i,j-1) + mu(i,j);
        line( [i,i], [l(i,j-1),l(i,j)],'color',col(j,:) )
    end
end
grid
box off

xlabel( 'data points i' )
ylabel( '\mu(i)')
leg = {}
for i = 1 : tsm.nc
    leg{i} = sprintf( '\\mu_{%d}', i );
end
title( 'Rule activation' )
legend( leg )
