%% Evaluate ARX model for given sched vars z == regressor x = u|y

% $Id$

function yp = tsm_evaluate_ARX( obj, u, y )
 
z = tsm_Regressor( u, obj.z_lag_u, y, obj.z_lag_y );
mu = obj.z_msf( z, obj.c, obj.m );
%Todo: nur wenn x_lag <> z_lag!
x = [z,ones(size(z,1),1)];

% Extendend Regressor
xe=[];
t = [1, size(x,2)];
for k=1:obj.nc
    xe = [xe repmat(mu(:,k), t) .* x];
end

yp = [ y(1:obj.maxlag); xe*obj.Theta ];

end
