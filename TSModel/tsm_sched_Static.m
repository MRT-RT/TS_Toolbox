%% Build the scheduling matrix z for static LiP models

%%
% Inputs:
%%
%  tsm      TS model
%  u        input matrix u(N x nu )
%  optional:
%  y        output vector y( N x 1 )
%
%%
% Outputs:
%%
%  z(N,nd)  scheduling matrix, nd=nu, ProductSpace: nd=nu+1
%
% $Id$

function z = tsm_sched_Static( tsm, u, y)

if nargin < 3
    y = tsm.y_ident;
    if nargin < 2
        u = tsm.u_ident;
    end
end

if tsm.ProductSpace
    z = [ u, y ];
else
    z = u;
end
