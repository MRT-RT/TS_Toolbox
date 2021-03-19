%% SSE: Compute Sum of Squared Errors (SSE)
%
% Inputs:
%   y_obsc observed values
%   y_pred predicted/model values
% Output: 
%   sse    Sum of Squared Errors (SSE)  

% $Id$

function sse = SSE( y_obsv, y_pred )

ec = ErrorCriteria( y_obsv, y_pred );

sse = ec.SSE;
