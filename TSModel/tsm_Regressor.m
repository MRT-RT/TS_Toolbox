%% tsm_Regressor: Build regressor matrix x for MISO data vectors (y|u)
%
% Inputs:
%  y(nr,1)            output vector
%  lag_y(1,nly)       lags in outputs 
%  optional:
%    u(nr,nu)         input vector
%    lag_u{nu}        lags in inputs { nu } 
%
% Outputs:
%  z(nr,lag_y+lag_u)  regressor matrix ( nr x lag_y,lag_u )
%
% $Id$


function x = tsm_Regressor(  u, lag_u, y, lag_y )

%<
%% Test data
if nargin == 0
    y = transpose( 1:10 );
    lag_y = [1,2,3];

    nu = 2;
    u = [-y , -y-0.2];
    lag_u = { [0], [0,2] };
end
%>

if ~isvector( y )
    error( 'tsm_Regressor: y not a vector')
end

if size(y,2) ~= 1
    warning( 'tsm_Regressor: y not a column vector')
    y = transpose(y);
end

ry = size( y, 1 ); % number of data points

if min( lag_y ) < 1
    warning( 'tsm_Regressor: lag(y) < 1')
end

maxlag =  max( lag_y );

%% optional input and lags given
if nargin > 2 
    
    [ru, nu] = size(u);
    %if nu ~= 1
    %    warning( 'tsm_Regressor: u not a column vector! MISO?')
    %end
    
    if nu == ry
        warning( 'tsm_Regressor: u not a n x  nu - matrix')
        u = transpose( u );
    elseif ru ~= ry
        error( 'tsm_Regressor: length u <> length y' )
    end
    
    if size( lag_u,2 ) ~= nu
        error( 'tsm_Regressor: cols lag_u <> nu' )
    end
    
    for iu = 1:nu
        maxlag = max( maxlag, max( lag_u{iu} ) );
    end
end

if maxlag > ry
    error( 'tsm_Regressor: max lag(u||y) > #rows(u/y) (not enough data)' )
end

% % for k=1:obj.length(fdata.n)
% %     XY(:,k)=y(maxlag-fdata.n(k)+1:end-fdata.n(k));
% % end
% % 
% % for k=1:length(fdata.m)
% %     XU(:,k)=u(maxlag-fdata.m(k)+1:end-fdata.m(k));
% % end

for k=1:length(lag_y)
    x(:,k) = y( maxlag-lag_y(k)+1 : end-lag_y(k) );
end

%% inputs given
if nu > 0
    for iu = 1 : nu
        lag = lag_u{iu};
        for k=1:length(lag)
            xu(:,k) = u( maxlag-lag(k)+1 : end-lag(k), iu );
        end
        x = [ x, xu ];
    end
end