%% Function to build the regressor nmatrix for LS models

% $Id$

function x = tsm_reg_ARX( tso, u, y )
if nargin > 1
    if size(u,2) ~= obj.nu
        error('tsm_reg_ARX: cols u <> n_u')
    end
else
     u = tso.u_ident;
end
if nargin > 2
    if size(y,2) ~= obj.ny
        error('tsm_reg_ARX: cols y <> n_y')
    end
else
    y = tso.y_ident;
end
if size(u,1) ~= size(y,1)
    error('tsm_reg_ARX: rows u <> rows y')
end

z = tsm_Regressor( u, tso.x_lag_u, y, tso.x_lag_y );

end
