function y = Friedman( x, n )

switch n
    case 2
        y = 10*sin( pi*x(:,1) .* x(:,2) );
    case 3
        y = 10*sin( pi*x(:,1).*x(:,2) ) + 20*( x(:,3)-0.5 ).^2;
    case 5
        y = 10*sin( pi*x(:,1).*x(:,2) ) + 20*( x(:,3)-0.5 ).^2 + 10*x(:,4) + 5*x(:,5);
end

