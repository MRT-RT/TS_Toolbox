%% Optimize LM parameters for ARX models
%
% $Id$

function y = tsm_optimize_ARX_MF( par, z, ts )

    % par0 = [ c | Theta ]

    t = ts.nA + ts.nB; 
    n = ts.nc * t;
    c = reshape( par(1:n), ts.nc, t );
    theta = par(n+1:end);

    mu = ts.z_msf( z, c, ts.m );
    x = [z, ones(size(z,1),1) ];
    xe=[];
    t = [1, size(x,2)];
    for k=1:ts.nc
        xe = [xe repmat(mu(:,k), t) .* x];
    end
    
    y = xe * transpose( theta );

end


