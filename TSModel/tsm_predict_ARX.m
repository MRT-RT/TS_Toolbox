%% Function tsm_predict_ARX: 
% Predict ARX model for given sched var z == regressor x = u|y
%%
% Inputs:
%   tso  TS model
%   u    input matrix u
%   y    output vector y 
% Output:
%   yp   predicted output

% $Id$

function yp = tsm_predict_ARX( tso, u, y )
 
z = tsm_Regressor( u, tso.z_lag_u, y, tso.z_lag_y );
mu = tso.z_msf( z, tso.v, tso.m );

%Todo: nur wenn x_lag <> z_lag!
x = [z,ones(size(z,1),1)];

% Extendend Regressor
xe=[];
t = [1, size(x,2)];
for k=1:tso.nv
    xe = [xe repmat(mu(:,k), t) .* x];
end

yp = [ y(1:tso.maxlag); xe*tso.Theta ];

end
