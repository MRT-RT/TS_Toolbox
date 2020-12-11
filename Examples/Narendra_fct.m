%% Narendra function

% $Id$

function [u,y] = Narendra_fct( n )

% $Id$

u = [];
while numel( u ) < n
    l = randi(20);
    a = 4*rand()-2;
    u = [u ; a*ones(l,1) ];
end
u = u(1:n);

y = zeros(n,1);
y(1) = 0;
y(2) = 0;
for k=2:n-1
    y(k+1) = ( y(k)*y(k-1)* (y(k)+2.5) ) / (1+y(k)^2+y(k-1)^2 ) +  u(k);
end