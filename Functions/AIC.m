%% AIC: Compute Akaike's Information Criteria (AIC)
%
% Inputs:
%   y_obsc observed values
%   y_pred predicted/model values
% Output:
%   aic    Akaike's Information Criteria (AIC)

% $Id$

function aic = AIC( y_obsv, y_pred, np )

ec = ErrorCriteria( y_obsv, y_pred, np );

aic = ec.AIC;
