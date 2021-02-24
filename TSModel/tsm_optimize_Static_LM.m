%% Optimize LM parameters for LiP models

% $Id$


function y = tsm_optimize_LiP_LM( par, u, ts )

    n1 = 1;
    n2 = ts.nc*ts.nu;
    B = reshape( par( n1:n2 ), ts.nc, ts.nu );
    n1 = n2+1;
    n2 = n1 + ts.nc-1;
    C = reshape( par( n1 : n2 ), ts.nc, 1 );

    mu = tsm_membership_FBF( u, ts.c, ts.m );
    ypl = transpose( ts.B * transpose(u) + ts.C );
    y = sum( mu .* ypl, 2 );

end


