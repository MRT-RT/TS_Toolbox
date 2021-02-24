%% Boxplot
function myBoxplot(x,xrange,y)

% $Id$

% width of box = 4 percent of x-range
dx = diff(xrange)/50;

qy = myQuantile(y,[0:0.25:1]);
%plot(x,qy(1),'r+',x,qy(5),'r+')
line([x,x]+[-dx,dx],[qy(1),qy(1)],'Color','k') % min
line([x,x],         [qy(1),qy(2)],'Color','k','LineStyle','--') % min->q1
line([x,x]+[-dx,dx],[qy(2),qy(2)],'Color','b') % q1
line([x,x]-dx,      [qy(2),qy(4)],'Color','b') % q1->q3
line([x,x]+dx,      [qy(2),qy(4)],'Color','b') % q1->q3
line([x,x]+[-dx,dx],[qy(3),qy(3)],'Color','r') % q2 = median
line([x,x]+[-dx,dx],[qy(4),qy(4)],'Color','b') % q3
line([x,x],         [qy(4),qy(5)],'Color','k','LineStyle','--') % q3->max
line([x,x]+[-dx,dx],[qy(5),qy(5)],'Color','k') % max
