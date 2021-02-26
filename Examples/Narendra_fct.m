%% Narendra function
%
% Create input $u$ as steps with width $l=1,\\ldots,20]$ for $N$ time steps 
% and compute output $y$ from the Narendra function 
%
% $Id$

function [ u, y ] = Narendra_fct( N )

% $Id$

u = [];
while numel( u ) < N
    l = randi(20);
    a = 4*rand()-2;
    u = [u ; a*ones(l,1) ];
end
u = u(1:N);

y = zeros(N,1);
y(1) = 0;
y(2) = 0;
for k=2:N-1
    y(k+1) = ( y(k)*y(k-1)* (y(k)+2.5) ) / (1+y(k)^2+y(k-1)^2 ) +  u(k);
end