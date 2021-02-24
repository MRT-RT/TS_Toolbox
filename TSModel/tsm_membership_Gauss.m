%% tsm_membership_Gauss
% Compute Gaussian membership degree mu( z )
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

[nc, nv ] = size(c );
n = size( z, 1 );

% Check nv  == size(z,2) == size(c,2)
if size(z,2) ~= nv 
    error( 'tsm_membership_Gauss: dim mismatch z / c' )
end

if size( z, 2 ) ~= nv
    error
mu = zeros(n,nv);
for ic = 1:nc
    dzc = sqrt( sum( (z - c(ic,:)).^2, 2 ) );
    mu(:,ic) = exp( m * dzc );
end

mu = mu ./ sum(mu,2);