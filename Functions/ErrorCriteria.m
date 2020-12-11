%% ErrorCtris: compute error criteria for two vectors
%
% Inputs:
%   y_t true values
%   y_p predicted/model values

% $Id$


function ec = ErrorCriteria( y_true, y_pred )

%% Check dimensions
[nrow,ncol] = size( y_true );
if min( [nrow,ncol] ) ~= 1
    error( 'ErrorCrits: y_true not n x1 vector')
end

if ~isequal( size(y_true), [nrow,ncol] )
    error( 'SimErrors: dim mismatch y_true <> y_pred')
end

dy = y_true - y_pred;

% Maximum absolute error (MAE)
ec.MAE = max( abs(dy) );

% Sum of squared errors (SSE)
ec.SSE = transpose( dy ) * dy;

% Mean squared error (MSE)
ec.MSE = ec.SSE / length(y_true);

% Root mean squared error (RMSE)
ec.RMSE = sqrt(ec.MSE);

% Normalized mean squared error (NMSE)
ym = mean( y_true );
dym = y_true - ym;
ec.NMSE = ec.MSE / ( transpose(dym) * dym);

% Best fit rate (BFR)
ec.BFR = 1 - sqrt(ec.NMSE);

