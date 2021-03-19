%% MSE: Compute Mean Squared Error (MSE)
%
% Inputs:
%   y_obsc observed values
%   y_pred predicted/model values
% Output:
%   mse    Mean Squared Error (MSE)

% Hint: ny = 1 due to MISO models in TS toolbox
% ts not given -> AIC/BIC = NaN

% $Id$

function mse = MSE( y_obsv, y_pred )

ec = ErrorCriteria( y_obsv, y_pred );

mse = ec.MSE;
