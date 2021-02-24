%% Optimize LM parameters for Static models

% $Id$


function y = tsm_optimize_Static_LM( par, u, ts )

    n1 = 1;
    n2 = ts.nv*ts.nu;
    B = reshape( par( n1:n2 ), ts.nv, ts.nu );
    n1 = n2+1;
    n2 = n1 + ts.nv-1;
    C = reshape( par( n1 : n2 ), ts.nv, 1 );

    mu = tsm_membership_FBF( u, ts.v, ts.m );
    ypl = transpose( ts.B * transpose(u) + ts.C );
    y = sum( mu .* ypl, 2 );

end


