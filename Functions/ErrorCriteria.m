%% ErrorCrtieria: compute error criteria for difference of two vectors
%
% Inputs:
%   y_obsc observed values
%   y_pred predicted/model values
% Optional:
%   np  number of parameters for computing AIC/BIC
% Outputs:
%   ec  struct with error criteria

% Hint: ny = 1 due to MISO models in TS toolbox
% ts not given -> AIC/BIC = NaN

% $Id$

function ec = ErrorCriteria( y_obsv, y_pred, np )

%% Check dimensions
if ~isvector( y_obsv )
    error( 'ErrorCrits: y_true not n x1 vector')
end

[nrow,ncol] = size( y_obsv );

if ~isequal( size(y_obsv), [nrow,ncol] )
    error( 'ErrorCriteria: dim mismatch y_obsv <> y_pred')
end

dy = y_obsv - y_pred;
N = length( dy );

% Maximum absolute error (MAE)
ec.MAE = max( abs(dy) );

% Sum of squared errors (SSE)
ec.SSE = transpose( dy ) * dy;

% Mean squared error (MSE)
ec.MSE = ec.SSE / N;

% Root mean squared error (RMSE)
ec.RMSE = sqrt( ec.MSE );

% Normalized mean squared error (NMSE)
ym = mean( y_obsv );
dym = y_obsv - ym;
ec.NMSE = ec.MSE / ( transpose(dym) * dym);

% Best fit rate (BFR)
ec.BFR = 1 - sqrt( ec.NMSE );

%% Akaiken and B criteria
if nargin > 2 && np > 0
    [ec.AIC,~,~,ec.BIC ] = aic( dym, N, np, 1 );
else
    ec.AIC = NaN;
    ec.BIC = NaN;
end
