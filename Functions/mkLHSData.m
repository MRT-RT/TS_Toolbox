%% mkLHSData
% $Id$

function lhsdata = mkLHSData(n,rows,cols,Plots)

if nargin < 1
    n=1
end

%V: c
if nargin < 2
    rows=50
end

% V : nD
if nargin < 3
    cols=2
end

if nargin < 4
    Plots = 0;
end

i=0;
lhsdata={};
if Plots
    c='krgbcmy'; lc=numel(c); ic=1;
    m='+x*sdo<>^>ph'; lm=numel(m); im=1;
    figure(1),clf
    hold all
end
l={};
for i=1:n
    lhsdata{i} = lhsdesign(rows,cols,'criterion','maximin') ;
    if Plots
    plot(lhsdata{i}(:,1),lhsdata{i}(:,2),[c(ic),m(im)])
    ic=ic+1;
    if ic>lc
        ic=1;
        im=im+1;
    end
    l{i} = sprintf('%3d',i);
    end
end
if Plots
    hold off
    legend(l,'Location','EastOutside')
end
save(sprintf('LHS_Data_rows%d_cols%d_n%d.mat',rows,cols,n),'lhsdata')
