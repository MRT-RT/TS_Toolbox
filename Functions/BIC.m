%% BIC: Compute Bayesian Information Criteria (BIC) 
%
% Inputs:
%   y_obsc observed values
%   y_pred predicted/model values
% Output:
%   aic    Bayesian Information Criteria (BIC)

% $Id$

function bic = BIC( y_obsv, y_pred, np )

ec = ErrorCriteria( y_obsv, y_pred, np );

bic = ec.BIC;
