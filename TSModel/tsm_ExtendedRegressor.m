%% Build extended regressor matrix phi = mu*[y|u|1] for ARX modells
%
%%
% Inputs:
%   mu(n,nx)      membership degrees of z
%   z(n,nx)       vector with output/input data
%%
% Outputs:
%   phi(n,nx+1)   extendend regressor matrix
%
% $Id$

function phi = tsm_ExtendedRegressor( mu, z )

%<
%% Test data
if nargin == 0
    n = 10;
    nx = 3;
    nc = 4;
    t=transpose(1:n);
    z = [ t, t+0.1, t+0.2 ];
    mu = rand(n,nx); 
    mu = mu ./ sum(mu,2);
    
end
%>

n = size( z,1 );         % number of data points
z = [ z, ones( n, 1 ) ]; % add column for affine term C
nx = size(z,2);          % number of regressor variables

nc = size(mu,2);
phi = zeros( n, nc*nx );

%% Build extended regressor matrix [ m_1*z | mu_2*z | ... | mu_nc*z ]
i1 = 1;
i2 = nx;
for k = 1:nc
    phi( : , i1:i2 ) =  repmat( mu(:,k),[1, nx] ).* z;
    i1 = i2 + 1;
    i2 = i2 + nx;
end
