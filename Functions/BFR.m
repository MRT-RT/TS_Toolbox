%% BFR: Compute Best Fit Rate (BFR)
%
% Inputs:
%   y_obsc observed values
%   y_pred predicted/model values
% Output:
%   bfr    Best fit rate (BFR)

% $Id$

function bfr = BFR( y_obsv, y_pred, np )

ec = ErrorCriteria( y_obsv, y_pred, np );

bfr = ec.BFR;
