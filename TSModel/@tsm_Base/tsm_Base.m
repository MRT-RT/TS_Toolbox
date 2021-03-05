%% TS-Toolbox: Definition of base class

% $Id$

classdef tsm_Base < handle & matlab.mixin.Copyable
    
    properties(Constant)
        Version = '0.91.'; % Version of TSModel class
    end
    
    properties
        Name = ''    % Name of object
        Date = ''    % Date of creation
        Comment = {} % Comments to object
        Debug = 0    % Debug level (0=none,1=flow,2=info)
    end
    
    methods
        
        function self = tsm_Base( name, comment )
            % constructor
            % input: name
            % input comment: initial comment, optional
            
            self.setName( name );
            self.Date = datetime( 'now' );
            if nargin > 1
                self.addComment( comment );
            end
        end
        
        function delete( self )
            % destructor
            self.debug( sprintf('tsm_Base/delete ''%s''', self.Name), 2 );
        end
        
        function setName( self, name )
            % set selfect name
            if ischar( name )
                self.Name = name;
            else
                error( 'tsm_Base/setName: argument not char' )
            end
        end
        
        function setDate( self, Date )
            if nargin < 2
                Date = datetime('now');
            end
            %todo: check if datestr
            self.Date = Date;
        end
        
        function self = addComment( self, comment )
            % add comment
            % input: comment: array co char
            if iscell( comment )
                for c = comment
                    self.Comment = [self.Comment, c];
                end
            else
                self.Comment = [self.Comment, comment];
            end
        end
        
        function setDebugLevel( self, level )
            % set debug level for diagnostic messages
            if isinteger( level ) and level < 0
                self.Debug = level;
            else
                error( 'tsm_Base/setDebug: argument not integer > 0' )
            end
        end

        function save( self, file )
            % save to MAT file
            % input: file - name of MAT file
            save( file, 'self' );
            self.debug( sprintf('tsm_Base: saved in ''%s''\n', file),2  );
        end
       
        function self = load( file )
            % load from MAT file
            % input: file - name of MAT file
            self = load( file, 'self'  );
            self.debug( sprintf('tsm_Base: loaded from ''%s''\n', file), 2 );
        end
        
        function disp( self )
            fprintf( 'Name: ''%s''\n', self.Name );
            fprintf( ' Type: ''%s''\n', class( self ) );
            if not( isempty( self.Date ) )
                fprintf( ' Date: ''%s''\n', self.Date );
            end
            if not( isempty( self.Comment ) )
                fprintf( ' Comments:\n' );
                for c = self.Comment
                    fprintf( '  ''%s''\n', c{1} );
                end
            end
        end
        
        function handle = plot( self )
            handle = None;
            warning( 'tsm_Base/plot: not implemented yet' )
        end
        
        function debug( self, message, level )
            % debug messages
            % input: level = integer > 0
            if level >= self.Debug
                fprintf( '%s\n', message )
            end
        end
        
    end
    
end
