%% Compute Fuzzy C-Means membership degree mu( z )
%
% $$ \mu_i(z) = \left[ \sum_{j=1}^c \left( \frac{||z-c_i||}{||z-c_j||}\right)^{m} \right]^{-1} $$
%
%%
% Inputs:
%%
%  z       scheduling vector (n x nv)
%  c       vector of clusters (nc xnv)
%  m       fuzziness parameter 2 / (nue-1)
%
%%
% Outputs:
%%
%  mu(n,nc) membership degree mu( z,c )

% $Id$

function mu = tsm_membership_FBF( z, c, m )

%%
n = size(z,1);
[nc,nv] = size(c);

%? Check nv  == size(z,2) == size(c,2)

%% dzc = || z - c(i,:) ||
dzc = zeros(n,nc);
for i=1:n
    for ic = 1:nc
        dzc(i,ic) = sqrt( sum( (z(i,:) - c(ic,:)).^2, 2 ) ).^m;
    end
end

%% sdzc = sum( 1/dzc^m )
sdzc = sum( 1./dzc, 2 );
% Test for reactivation --> sum == 0
tsdzc = sdzc < realmin;
if any( tsdzc )
    warning( 'membership degree: reactivation?' )
end

%% mu( z )
mu = zeros(n,nc);
for ic=1:nc
    mu(:,ic) = 1 ./ ( dzc(:,ic) .* sdzc );
    % z = c(i)
    mu( dzc(:,ic) < realmin ,ic ) = 1;
    % reactivation
    mu( tsdzc, ic ) = 1/nc;
end
