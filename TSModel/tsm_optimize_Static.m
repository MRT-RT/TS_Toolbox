%% Optimize MF&LM parameters for Static models

% $Id$


function y = tsm_optimize_Static( par, u, ts )

% par0 = [ reshape(obj.c,1,obj.nv*obj.nu),...
%          reshape(obj.B,1,obj.nv*obj.nu),...
%          reshape(obj.C,1,obj.nv) ];

    n2 = ts.nv*ts.nu;
    v = reshape( par(1:n2), ts.nv, ts.nu );
    n1 = n2+1;
    n2 = n2 + n2;
    B = reshape( par( n1:n2 ), ts.nv, ts.nu );
    n1 = n2+1;
    n2 = n1 + ts.nv-1;
    C = reshape( par( n1 : n2 ), ts.nv, 1 );

    mu = ts.z_msf( u, v, ts.m ); 
    ypl = transpose( B * transpose(u) + C );
    y = sum( mu .* ypl, 2 );

end


