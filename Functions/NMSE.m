%% NMSE: Compute Normalized Mean Squared Error (NMSE)
%
% Inputs:
%   y_obsc observed values
%   y_pred predicted/model values
% Output:
%   nmse   Normalized Root Mean Squared Error (NMSE)

% $Id$

function nmse = NMSE( y_obsv, y_pred )

ec = ErrorCriteria( y_obsv, y_pred );

nmse = ec.NMSE;
