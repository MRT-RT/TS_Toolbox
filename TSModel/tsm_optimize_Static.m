%% Optimize MF&LM parameters for LiP models

% $Id$


function y = tsm_optimize_LiP( par, u, ts )

% par0 = [ reshape(obj.c,1,obj.nc*obj.nu),...
%          reshape(obj.B,1,obj.nc*obj.nu),...
%          reshape(obj.C,1,obj.nc) ];

    n2 = ts.nc*ts.nu;
    c = reshape( par(1:n2), ts.nc, ts.nu );
    n1 = n2+1;
    n2 = n2 + n2;
    B = reshape( par( n1:n2 ), ts.nc, ts.nu );
    n1 = n2+1;
    n2 = n1 + ts.nc-1;
    C = reshape( par( n1 : n2 ), ts.nc, 1 );

    mu = tsm_membership_FBF( u,c,ts.m );
    ypl = transpose( B * transpose(u) + C );
    y = sum( mu .* ypl, 2 );

end


