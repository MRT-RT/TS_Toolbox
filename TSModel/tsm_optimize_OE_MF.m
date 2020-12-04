%% Optimize MF&LM parameters for LIP models
%
% $Id$


function y = tsm_optimize_OE_MF( par, ts )

% par0 = [ c | Theta ]

    t = ts.nA + ts.nB; 
    n = ts.nc * t;
    c = reshape( par(1:n), ts.nc, t );
    theta = transpose( par(n+1:end) );

    y = tsm_evaluate_OE( ts, ts.u_ident, ts.y_ident, theta ) - ts.y_ident;

end


