%% tsm_predict_oe
% Predict OE model for given sched vars z == regressor x = u|y
%%
%  Inputs:
%   tso   TS model      (instance of class TSModel)
%   u     input vector  (n x nu )
%   y     initial output vector  (maxlag x 1 )
%   Theta parameter od local modells
%  Outputs:
%   yp     predicted output vector (n x 1 )

% $Id$

function yp = tsm_predict_OE( tso, u, y, Theta )

if nargin < 4
    Theta = tso.Theta;
end

n = size(y,1);
yp = zeros( n, 1);
yp(1:tso.x_maxlag) = y(1:tso.x_maxlag);

ny = length(tso.x_lag_y);
x = zeros( 1, tso.nx );
for i=tso.x_maxlag+1:n

    for k=1:ny
        x(k) = yp( i+1-tso.x_lag_y(k) );
    end
    
    for iu=1:tso.nu
        lag = tso.x_lag_u{iu};
        nu = length( lag );
        for k = 1:nu
            x(ny+k) = u( i-lag(k) );
        end
    end
    mu = tso.z_msf( x, tso.v, tso.m );
    xe=[];
    t = [1, size(x,2)+1];
    for k = 1:tso.nv
        xe = [xe repmat(mu(:,k), t) .* [x,1]];
    end
    
    yp(i) = xe*Theta;
    
end

