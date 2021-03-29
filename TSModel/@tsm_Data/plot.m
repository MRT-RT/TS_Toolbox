%% Plot tsm_Data object
%
% $Id$

function handle = plot( self, fig )

if nargin > 1
    handle = figure( fig );
else
    handle = figure;
end

for i=1:self.nu
    subplot( self.nu + 1, 1, i );
    plot( self.t,self.u(:,i) )
    if i == 1
        title( sprintf( 'Dataset: ''%s'' data-points:%d', self.Name, self.n ) )
    end
    if isempty( self.u_limit ) 
        ylim( [min(self.u(:,i)),max(self.u(:,i))] )
    else
        ylim( self.u_limit(i,:) )
    end
    if isempty( self.t_limit )
        xlim( [self.t(1),self.t(end)] ),
    else
        xlim( self.t_limit ),
    end
    grid on
    
    set( gca, 'xticklabel',{} )
    ylabel( self.u_label{i} )
end
subplot( self.nu + 1, 1, self.nu + 1 );
plot( self.t, self.y  )
if isempty( self.y_limit )
    ylim( [min(self.y),max(self.y)] )
else
    ylim( self.y_limit )
end
if isempty( self.t_limit )
    xlim( [self.t(1),self.t(end)] ),
else
    xlim( self.t_limit ), 
end
grid
ylabel( self.y_label )
xlabel( self.t_label )
%? Limits