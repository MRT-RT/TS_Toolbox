%% function bic: Bayesian Information Criterion

% $Id$

function BIC = bic( E, N, np, ny )

%     'BIC': Bayesian Information Criterion
%            BIC = N*log(det(E'*E/N)) + N*(ny*log(2*pi) + 1) + d*log(N)

[ ~,~,~,BIC ] = aic( E, N, np, ny );

