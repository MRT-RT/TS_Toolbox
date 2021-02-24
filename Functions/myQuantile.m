%% myQuantile
%%
% Inputs
% v vector
% p quantile percentages (optional -> q = 0/0.25(0.5/0.75/1)
%%
% Outputs
% q quantile(p)
%
% $Id$

function q = myQuantile( v, p )

% p=0:    Q0 min
% p=0.25: Q1
% p=0.5:  Q2 median
% p=0.75: Q3 
% p=1:    Qn max

% all quantiles
if nargin < 2
    p = 0:0.25:1;
end

if size(v,1) < size(v,2)
    v = transpose(v);
end
nv = size(v,1);
v = sort(v);

if size(p,1) < size(p,2)
    p = transpose(p);
end
np = size(p,1);
q = zeros(np,1);

ind = floor( nv*p );
for ip=1:np
    if ind(ip) <= 1
        % minimal value
        q(ip) = v(1);
    elseif ind(ip) >= nv
        % maximal value
        q(ip) = v(end);
    elseif mod( ind(ip),nv )
        q(ip) = (v(ind(ip)) + v(ind(ip)+1) ) /  2;
    else
        q(ip) = v(ip);
    end
end

