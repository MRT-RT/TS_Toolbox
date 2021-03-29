%% Display tsm_Data object
%
% $Id$

function disp( self )

disp@tsm_Base( self )
fprintf( '#data-points n = %d\n', self.n );
if self.dt > 0
    fprintf( 'time dt = %g', self.dt );
    if self.t_limit
        fprintf( ' limit: [%g, %g]', self.t_limit );
    end
    if self.t_unit
        fprintf( ' unit: ''%s''', self.t_unit );
    end
    if self.t_label
        fprintf( ' label: ''%s''', self.t_label );
    end
    fprintf( '\n' )
end
fprintf( '#inputs u = %d\n', self.nu );
for i=1:self.nu
    fprintf( ' u%d: ',i )
    if not(isempty( self.u_limit ) )
        fprintf( ' limit: [%g, %g]', self.u_limit(i,:) );
    end
    if not(isempty( self.u_unit ) )
        fprintf( ' unit: ''%s''', self.u_unit{i} );
    end
    if not(isempty( self.u_label ) )
        fprintf( ' label: ''%s''', self.u_label{i} );
    end
    fprintf( '\n' )
end
fprintf( '#outputs y = %d', self.ny );
if self.y_limit
    fprintf( ' limit: [%g, %g]', self.y_limit );
end
if self.y_unit
    fprintf( ' unit: ''%s''', self.y_unit );
end
if self.y_label
    fprintf( ' label: ''%s''', self.y_label );
end
fprintf( '\n' )
