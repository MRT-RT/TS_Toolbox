%% Build the scheduling matrix z for ARX models

% $Id$

function z = tsm_sched_ARX( tso, u, y)
if nargin > 1
    if size(u,2) ~= tso.nu
        error('tsm_Sched_ARX: cols u <> n_u')
    end
else
    u = tso.u_ident;
end
if nargin > 2
    if size(y,2) ~= tso.ny
        error('tsm_Sched_ARX: cols y <> n_y')
    end
else
    y = tso.y_ident;
end
if size(u,1) ~= size(y,1)
    error('tsm_Sched_ARX: rows u <> rows y')
end

z = tsm_Regressor( u, tso.z_lag_u, y, tso.z_lag_y );

end
