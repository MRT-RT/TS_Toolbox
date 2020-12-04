%% Evaluate LiP model for given sched vars z == regressor x

function [ yp, ypl ] = tsm_evaluate_LiP( obj, z )
mu = obj.z_msf( z, obj.c, obj.m );
ypl = transpose( obj.B * transpose(z) + obj.C );
yp = sum( mu .* ypl, 2 );
end
