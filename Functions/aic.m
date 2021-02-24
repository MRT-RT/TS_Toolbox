%% function aic: Akaike's Information Criterion (AIC)

% $Id$

function [AIC, nAIC, AICc, BIC] = aic( E, N, np, ny )

% 'nAIC': Normalized Akaike's Information Criterion
%            nAIC = log(det(E'*E/N)) + 2d/N
% 'AIC': Akaike's Information Criterion (raw measure)  (default)
%            AIC = N*nAIC + N*(ny*log(2*pi) + 1)
% 'AICc': AIC corrected for small sample size
%            AICc = AIC + d*(d+1)/(N-d-1))
% 'BIC': Bayesian Information Criterion
%            BIC = N*log(det(E'*E/N)) + N*(ny*log(2*pi) + 1) + d*log(N)

if isvector(E)
    [r,c] = size( E );
    if r == 1 
        warning( 'aic: arg E not a column vector->transpose(E)' )
        E = transpose( E );
        NE = c;
    else
        NE =r;
    end
else
    error( 'aic: arg E not a vector!' )
end

if NE <= 1 || NE > N
  warning( 'aic: length(E) not 1 ... N' )
end  

nAIC = log( det( E'*E/N ) ) + 2*np/N;
AIC = N * (nAIC + (ny*log(2*pi) + 1) );
AICc = AIC + np*(np+1)/(N-np-1);
BIC = AIC + N*(ny*log(2*pi) + 1) + np*log(N);

