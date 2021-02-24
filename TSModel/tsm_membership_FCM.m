%% tsm_membership_FCM
% Compute Fuzzy C-Means membership degree mu( z )
%
% $$ \mu_i(z) = \left[ \sum_{j=1}^n_v \left( \frac{||z-v_i||}{||z-v_j||}\right)^{m} \right]^{-1} $$
%
% Set mu = 1, if ||z_j-v_i|| < realmin
% Set mu = 0, if sum ||z_j-v_i|| < realmin
%
%%
% Inputs:
%%
%  z        scheduling vector (N x nd)
%  v        vector of clusters (nv x nd)
%  m        fuzziness parameter m = 2 / (nue-1)
%
% Outputs:
%%
%  mu(N,nd)  membership degree mu( z,v )
%
% $Id$

function mu = tsm_membership_FCM( z, v, m )

N = size( z, 1 );
[ nv, nd ] = size( v );

% Check nv  == size(z,2) == size(v,2)
if size( z, 2 ) ~= nd 
    error( 'tsm_membership_FBF: dim mismatch z / v' )
end

% Distance of all data z to cluster centers c: dzc = || z - v(i,:) ||
dzv = zeros(N,nv);
for i=1:N
    for iv = 1:nv
        dzv(i,iv) = sqrt( sum( (z(i,1:nd) - v(iv,:)).^2, 2 ) ).^m;
    end
end

% Inverse sum of distances sdzv = sum( 1/dzv^m )
sdzv = sum( 1./dzv, 2 );

% Search for reactivation --> sum ~ 0
tsdzv = sdzv < realmin;
if any( tsdzv )
    warning( 'tsm_membership_FBF: possible reactivation' )
end

% membership degree mu( z )
mu = zeros( N, nv );
for iv=1:nv
    mu(:,iv) = 1 ./ ( dzv(:,iv) .* sdzv );
    % z = c(i)
    mu( dzv(:,iv) < realmin ,iv ) = 1;
    % reactivation
    mu( tsdzv, iv ) = 1/nv;
end
