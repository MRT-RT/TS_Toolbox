%% Evaluate ARX model for given sched vars z == regressor x = u|y

%%
% Inputs:
%   obj:  TS model      (instance of class TSModel)
%   u:    input vector  (n x 1 )
%   y:    output vector  (n x 1 )
%%
% Outputs:
%   yp prdicted output vector (n x 1 )

% $Id$



function yp = tsm_evaluate_OE( obj, u, y, Theta )

if nargin < 4
    Theta = obj.Theta;
end

n = size(y,1);
yp = zeros( n, 1);
yp(1:obj.x_maxlag) = y(1:obj.x_maxlag);

ny = length(obj.x_lag_y);
x = zeros( 1, obj.nx );
for i=obj.x_maxlag+1:n

    for k=1:ny
        x(k) = yp( i+1-obj.x_lag_y(k) );
    end
    
    for iu=1:obj.nu
        lag = obj.x_lag_u{iu};
        nu = length( lag );
        for k = 1:nu
            x(ny+k) = u( i-lag(k) );
        end
    end
    mu = obj.z_msf( x, obj.c, obj.m );
    xe=[];
    t = [1, size(x,2)+1];
    for k = 1:obj.nc
        xe = [xe repmat(mu(:,k), t) .* [x,1]];
    end
    
    yp(i) = xe*Theta;
    
end

