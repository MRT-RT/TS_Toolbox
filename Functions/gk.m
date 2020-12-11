%% Clustering with fuzzy covariance matrix (Gustafson-Kessel algorithm)
%
%--------------------------------------------------
% Inputs:
% z ... N by n data matrix
% nc ... number of clusters
% m ... fuzziness exponent (m > 1)
% tol ... termination tolerance (tol > 0)
%--------------------------------------------------
% Outputs:
% c ... cluster means (centers)
% U ... fuzzy partition matrix
% F ... cluster covariance matrices
%----------------- prepare matrices ----------------------------------

% $Id$

function [c,U,F] = gk(z,nc,m,tol)

[N,n] = size(z);   % data size
N1 = ones(N,1); 
n1 = ones(n,1); 
c1 = ones(1,nc);   % aux. variables
U = zeros(N,nc);   % partition matrix
d = U;             % distance matrix
F = zeros(n,n,nc); % covariance matrix
%----------------- initialize U --------------------------------------
minZ = c1'*min(z); 
maxZ = c1'*max(z); % data limits
c = minZ + (maxZ - minZ).*rand(nc,n); % random centers
for j = 1 : nc
    Zc = z - N1*c(j,:);
    d(:,j) = sum((Zc.^2),2);   % distances
end
d = (d+eps).^(-1/(m-1));       % inverse dist.
U0 = (d ./ (sum(d,2)*c1));     % part. matrix
%----------------- iterate --------------------------------------------
while max( max(U0-U) ) > tol   % no convergence
    U = U0; 
    Um = U.^m; sumU = sum(Um); % aux. vars
    c = (Um'*z)./(n1*sumU)';   % clust. centers
    for j = 1 : nc             % for all clusters
        Zc = z - N1*c(j,:);    % auxiliary var
        f = n1*Um(:,j)'.*Zc'*Zc/sumU(j); % cov. matrix
        d(:,j) = sum( Zc*( det(f)^(1/n)*inv(f) ).*Zc,2 ); % distances
    end
    % d = (d+1e-100).^(-1/(m-1));   % inverse dist.
    d = (d+eps).^(-1/(m-1));      % inverse dist.
    U0 = (d ./ (sum(d,2)*c1));    % part. matrix
end
%----------------- create final F and U -------------------------------
U = U0; 
Um = U.^m; 
sumU = n1*sum(Um);
for j = 1 : nc
    Zc = z - N1*c(j,:);
    F(:,:,j) = n1*Um(:,j)'.*Zc'*Zc/sumU(1,j);
end
%----------------- end of function ------------------------------------
