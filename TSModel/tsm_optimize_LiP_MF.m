%% Optimize MF parameters for LiP models

% $Id$

function y = tsm_optimize_LiP_MF( par, u, ts )

    c = reshape( par(1:ts.nc*ts.nu), ts.nc, ts.nu );

    mu = tsm_membership_FBF( u,c,ts.m );
    ypl = transpose( ts.B * transpose(u) + ts.C );
    y = sum( mu .* ypl, 2 );

end


