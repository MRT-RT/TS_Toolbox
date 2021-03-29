%% fcm_Euclidean
% The function fcm_Euclidean performs Fuzzy-C-Means clustering with an
% euclidean distance measure.

% $Id$

function [v, distout, J ] = fcm_Euclidean(data, cluster, options)


% get parameters
if nargin > 2
    m = options(1);
    max_iter = options(2);
    min_impro = options(3);
    display = options(4);
else % no options given
    m = 2;
    max_iter = 100;
    min_impro = 1e-5;		% Min. improvement
    display = 1;           	% info display during iteration
end

[N,n] = size(data);
[Nc,nc] = size(cluster);
X1 = ones(N,1);

if max(Nc,nc) == 1 		% only number of cluster given
    
    c = cluster;
    mm = mean(data);             %mean of the data (1,n)
    aa = max(abs(data - ones(N,1)*mm)); %
    v = 2*(ones(c,1)*aa).*(rand(c,n)-0.5) + ones(c,1)*mm;
    for j = 1 : c
        xv = data - X1*v(j,:);
        d(:,j) = sum((xv*eye(n).*xv),2);
    end
    d = (d+1e-10).^(-1/(m-1));
    cluster = (d ./ (sum(d,2)*ones(1,c)));
    
else
    
    c = size(cluster,2);
    fm = cluster.^m; sumf = sum(fm);
    v = (fm'*data)./(sumf'*ones(1,n)); %
    
end

f = zeros(N,c);                 % partition matrix

% Iterate
iter = 0;                       % iteration counter
while iter <= max_iter
    
    iter = iter + 1;
    
    if abs( max(max(cluster - f)) ) < min_impro
        break
    end
    
    f = cluster;
    % Calculate centers
    fm = f.^m;
    sumf = sum(fm);
    v = (fm'*data) ./ (sumf'*ones(1,n));
    for j = 1 : c
        xv = data - X1*v(j,:);
        d(:,j) = sum((xv*eye(n).*xv),2);
    end
    distout=sqrt(d);
    J(iter) = sum(sum(cluster.*d));
    
    if display
		fprintf('Iteration count = %d, obj. fcn = %f\n', iter, J(iter));
    end
    
    % Update f0
    d = (d+1e-10).^(-1/(m-1));
    cluster = (d ./ (sum(d,2)*ones(1,c)));
end

