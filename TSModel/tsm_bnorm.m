%% Norm of vector = z^T * z

function n = tsm_bnorm( z )

if size(z,1) > size(z,2)
    n = z'*z;
else
    n = z*z';
end

