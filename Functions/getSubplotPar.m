%% Compute number of rows/colums for subplots
%
% $Id$

function [rows,cols] = getSubplotPar( n, dim )

if nargin < 2
    dim = 2;
end

n = factorial( n  ) / (factorial( n-dim ) * dim );
rows = round( sqrt( n ) );
cols = ceil( n / rows );
