%% MAE: Compute Maximum Absolute Error (MAE)
%
% Inputs:
%   y_obsc observed values
%   y_pred predicted/model values
% Output:
%   mae    Maximum Absolute Error (MAE)

% $Id$

function mae = MAE( y_obsv, y_pred )

ec = ErrorCriteria( y_obsv, y_pred );

mae = ec.MAE;
