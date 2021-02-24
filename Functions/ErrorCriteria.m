%% ErrorCrtieria: compute error criteria for difference of two vectors
%
% Inputs:
%   y_t true values
%   y_p predicted/model values
% Optional:
%   np  number of parameters AIC/BIC
% Outputs:
%   ec  struct with error criteria

% Hint: ny = 1 due to MISO models in TS toolbox
% ts not given -> AIC/BIC = -1

% $Id$

function ec = ErrorCriteria( y_true, y_pred, np )

%% Check dimensions
if ~isvector( y_true )
    error( 'ErrorCrits: y_true not n x1 vector')
end

[nrow,ncol] = size( y_true );

if ~isequal( size(y_true), [nrow,ncol] )
    error( 'ErrorCrits: dim mismatch y_true <> y_pred')
end

dy = y_true - y_pred;
N = length( dy );

% Maximum absolute error (MAE)
ec.MAE = max( abs(dy) );

% Sum of squared errors (SSE)
ec.SSE = transpose( dy ) * dy;

% Mean squared error (MSE)
ec.MSE = ec.SSE / N;

% Root mean squared error (RMSE)
ec.RMSE = sqrt(ec.MSE);

% Normalized mean squared error (NMSE)
ym = mean( y_true );
dym = y_true - ym;
ec.NMSE = ec.MSE / ( transpose(dym) * dym);

% Best fit rate
ec.BFR = 1 - sqrt(ec.NMSE);

%% TS model for numer of parameters
if np > 0
    [ec.AIC,~,~,ec.BIC ] = aic( dym, N, np, 1 );
else
    ec.AIC = -1;
    ec.BIC = -1;
end
