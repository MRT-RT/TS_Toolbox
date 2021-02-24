function h = myVoronoi( x,c,l,w )

% $Id$

if nargin < 3
    l = '--';
end
if nargin < 2
    c = 'r';
end
if nargin < 4
    w = 1;
end
[xv,yv]=voronoi(x(:,1),x(:,2));
h = plot(xv,yv,'Color',c,'LineStyle',l,'LineWidth',w);
