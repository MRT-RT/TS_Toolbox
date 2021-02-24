%% Optimize MF & LM parameters for ARX models

% $Id$

function y = tsm_optimize_ARX( par, z, ts )

    % par0 = [ v | Theta ]
    
    %ad: size(v) = nv x (z_lag_u + z_lag_y)
    %ad: size(Theta) = nv x (nA + nB + 1)
    
    t = ts.nA + ts.nB; 
    n = ts.nv * t;
    v = reshape( par(1:n), ts.nv, t );
    theta = par(n+1:end);

    mu = ts.z_msf( z, v, ts.m );
    x = [z, ones(size(z,1),1) ];
    xe=[];
    t = [1, size(x,2)];
    for k=1:ts.nv
        xe = [xe repmat(mu(:,k), t) .* x];
    end
    
    y = xe * transpose( theta );

end


