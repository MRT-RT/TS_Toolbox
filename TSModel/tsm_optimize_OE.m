%% Optimize MF&LM parameters for LIP models

% $Id$


function y = tsm_optimize_OE( par, ts )

% par0 = [ v | Theta ]

    t = ts.nA + ts.nB; 
    n = ts.nv * t;
    v = reshape( par(1:n), ts.nv, t );
    theta = transpose( par(n+1:end) );

    y = tsm_predict_OE( ts, ts.u_ident, ts.y_ident, theta ) - ts.y_ident;

end


