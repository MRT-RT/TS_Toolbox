%% makeGrid
%
% Generate combinations n x nDim grid with nx points per dimension => nx^nDim points
%%
% Inputs:
%  n       number of data points per 
%  nz      number of scheduling variables 
%  limits  
%%
% Outputs:
%  x
%
% $Id$

function [ Xg, nXg ] = makeGrid( n, nz, limits )

%nxg = floor( nthroot(nX,nDim) );

Xg = zeros(n^nz,nz);
s1 = '';
s2 = ' ';
for i = 1:nz
    s1 = sprintf('%s %clinspace(%g,%g,%d)',...
        s1,s2,limits(i,1),limits(i,2),n );
    s2=',';
end
eval( sprintf(' Xg = transpose( combvec(  %s ) );',s1 ) )
nXg = size( Xg,1);
