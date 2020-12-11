%% Norm of vector = z^T * z

% $Id$

function n = bnorm( z )

if size(z,1) > size(z,2)
    n = z'*z;
else
    n = z*z';
end

