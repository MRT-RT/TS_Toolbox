%% RMSE: Compute Root Mean Squared Error (RMSE)
%
% Inputs:
%   y_obsc observed values
%   y_pred predicted/model values
% Output:
%   rmse   Root Mean Squared Error (RMSE)

% $Id$

function rmse = RMSE( y_obsv, y_pred )

ec = ErrorCriteria( y_obsv, y_pred );

rmse = ec.RMSE;
