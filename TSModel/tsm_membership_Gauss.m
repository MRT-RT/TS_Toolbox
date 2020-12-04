%% Compute Gaussian membership degree mu( z )
%%
% 
% $$\mu_i(z) = e^{m\cdot ||z-c_i||}$$
 

%%
% Inputs:
%%
%   z    scheduling variables matrix (n x nz )
%   c    cluster centers (nc x nz )
%   m    fuzziness exponent (m > 0)
%%
% Outputs:
%%
%   mu   membership degrees mu(z) (n x nc)
%%
% $Id$

function mu = tsm_membership_Gauss( z, c, m)

% Test data
if nargin == 0
    m = 2 / (1.2-1);
    n = 10;
    nc = 3;
    nv = 3;
    z = zeros( n, nv );
    t  = transpose( 1: n );
    for i=1:nv
        z (:,i) = t + i/10;
    end
    x = randi( n,nc,1 );
    c = z( x, : );
end

[n,nv] = size( z );
% nv == size(c,2) 

mu = zeros(n,nv);
for ic = 1:nc
    dzc = sqrt( sum( (z - c(ic,:)).^2, 2 ) );
    mu(:,ic) = exp( m * dzc );
end

mu = mu ./ sum(mu,2);