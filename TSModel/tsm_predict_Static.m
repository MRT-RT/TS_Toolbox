%% Pretict static model for given sched vars z == regressor x

%%
% Inputs: 
% tsm   TS model
% z     scheduling matrix
%%
% Outputs:
% yp    predicted output
% optional
% ypl   predicted outputs of local models

% $Id$

function [ yp, ypl ] = tsm_predict_Static( tsm, z )

    mu = tsm.z_msf( z, tsm.v, tsm.m );
    ypl = transpose( tsm.B * transpose(z(:,1:tsm.nu)) + tsm.C );
    yp = sum( mu .* ypl, 2 );

end
