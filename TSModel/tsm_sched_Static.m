%% Build the scheduling matrix z for linear static (LS) models

% $Id$

function z = tsm_sched_LiP( tso, u, y)

if nargin == 2
    z = u;
elseif nargin == 3
    z = [u,y];
else
    if tso.ProductSpace
        z = [tso.u_ident,tso.y_ident];
    else
        z = tso.u_ident;
    end
end
