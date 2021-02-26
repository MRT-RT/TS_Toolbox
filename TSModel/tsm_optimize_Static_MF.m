%% Optimize MF parameters for Static models

% $Id$

function y = tsm_optimize_Static_MF( par, u, ts )

    v = reshape( par(1:ts.nv*ts.nu), ts.nv, ts.nu );

    mu = ts.z_msf( u, v, ts.m );     
    ypl = transpose( ts.B * transpose(u) + ts.C );
    y = sum( mu .* ypl, 2 );

end


