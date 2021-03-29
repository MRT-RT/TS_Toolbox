%% Class definition TS-Toolbox: input/output data

% $Id$
%{
    MISO dataset
    t(n x 1)  (optional dt)
    u( n x nu)
    y ( n x 1)

    ToDo:
        split dataset for ident/validate
%}

classdef tsm_Data < tsm_Base
    
    properties
        
        n = 0   % number rof rows
        nu  = 0 % number of inputs (=columns in u)
        ny  = 1 % number of outputs (=columns in y)
        t = []  % time vector
        u = []  % array of input vectors n x nu
        y = []  % output vector n x 1
        
        dt = -1 % time delta
        t_label = 't'
        u_label = {}
        y_label = 'y'
        
        t_unit = '' % Unit of time vector
        u_unit = {} % Units of input vectors
        y_unit = '' % Unit of output vector
        
        t_limit = [] % Limits of time axis t_limit(1,2)
        u_limit = [] % Limits of input axss u_limit(nu,2)
        y_limit = [] % Limits of output axis y_limit(1,2)
        
    end
    
    methods
        
        function self = tsm_Data( name, u, y, varargin )
            
            p = inputParser;
            valCoC = @(x) iscell(x) || ischar( x );
            p.addRequired( 'name', @ischar )
            p.addRequired( 'u', @ismatrix )
            p.addRequired( 'y', @isvector )
            p.addParameter('t',[],@isvector)
            p.addParameter('comment',{},valCoC)
            p.addParameter('dt',-1,@isscalar)
            p.addParameter('tlabel','',@ischar)
            p.addParameter('ulabel',{},@iscell)
            p.addParameter('ylabel','',@ischar)
            p.addParameter('tunit','',@ischar)
            p.addParameter('uunit',{},@iscell)
            p.addParameter('yunit','',@ischar)
            p.addParameter('tlim', [],@isvector)
            p.addParameter('ulim', [],@ismatrix)
            p.addParameter('ylim', [],@isvector)
            p.parse( name, u, y, varargin{:} )
            opts = p.Results;
            
            self@tsm_Base( name );
            
            self.set_u( u );
            self.set_y( y )
            if not( isempty( opts.comment ) )
                self.addComment( opts.comment )
            end
            
            if ~isempty( opts.t )
                set_t( opts.t );
            end
            if opts.dt > 0
                self.dt = opts.dt;
            end
            
            if isempty( opts.tlabel )
                self.t_label = 't';
            else
                self.set_Label( 't',opts.tlabel )
            end
            if isempty( opts.tlabel )
                for i=1:self.nu
                    self.u_label{i} = sprintf( 'u%d',i);
                end
            else
                self.set_Label( 'u',opts.ulabel )
            end
            if isempty( opts.ylabel )
                self.t_label = 'y';
            else
                self.set_Label( 'y',opts.ylabel )
            end
            
            self.set_Unit( 't',opts.tunit )
            self.set_Unit( 'u',opts.uunit )
            self.set_Unit( 'y',opts.yunit )
        end
        
        function v = make_ColVector( self, v )
            [ r, c ] =  size( v );
            if r == 1 % a row vector
                v = transpose( v );
                n = c;
            elseif c == 1 % not a col vector
                n = r;
            else
                error( 'tsm_Data/make_ColVector: v not a vector (%d x %d)',r,c )
            end
            if self.n > 0 && not( n == self.n )
                error( 'tsm_Data/make_ColVector: len( v ) inconsitent (%d : %d)',n,self.n )
            end
        end
        
        % set time vector t( n x 1 )
        function set_t( self, t )
            
            t = self.make_ColVector( t );
            
            % Check consistent sample time
            dt = t(2) - t(1);
            if  self.dt > 0 && dt ~= self.dt
                error( 'tsm_Data: t (%g) and different dt (%g) given', ...
                    dt,self.dt )
            else
                self.dt = dt;
                self.t = t;
            end
            
        end
        
        % set input vectors u( n x nu )
        function  set_u( self, u)
            if ismatrix( u )
                [ n, nu ] = size( u );
                % nr=n or nc=n ???
                if self.n > 0 && self.n ~= n
                    error( 'tsm_Data/set_u: len u(%d) not n=%d)', n,self.n )
                else
                    self.n = n;
                end
                self.nu = nu;
                self.u = u;
            else
                self.nu = 0;
                self.u = [];
            end
        end
        
        % output vector
        function set_y( self, y )
            y = self.make_ColVector( y );
            n = size( y, 1 );
            if self.n > 0 && n ~= self.n
                error( 'tsm_Data/set_y: len y(%d) not n=%d', n,self.n )
            else
                self.n = n;
            end
            self.y = y;
        end
        
        % Label for
        function set_Label( self, v1, v2 )
            
            if not( isempty( v1 ) )
                
                if nargin == 2 % all 3 labels given { t,{u},y }
                    
                    if ~iscell( v1 )
                        error( 'tsm_Data.set_Label arg 1 not cell' )
                    end
                    if length( v1 ) ~= 3
                        error( 'tsm_Data.set_Label len(arg 1) <> 3' )
                    end
                    if not( ischar( v1{1} ) )
                        error( 'tsm_Data.set_Label: label t not string' )
                    end
                    if length( v1{3} ) ~= self.nu
                        error( 'tsm_Data.set_Label: label u not string(%d)',self.nu )
                    end
                    for i=1:self.nu
                        if not( ischar( v1{2}{i} ) )
                            error( 'tsm_Data.set_Label: label u(%d) not string',i )
                        end
                    end
                    if not( ischar( v1{3} ) )
                        error( 'tsm_Data.set_Label: label y not string' )
                    end
                    self.t_label = v1{1};
                    self.u_label = v1{2};
                    self.y_label = v1{3};
                    
                elseif nargin == 3 % var, label
                    
                    if v1 == 't'
                        if not( ischar( v2 ) )
                            error( 'tsm_Data.set_Label t not string' )
                        else
                            self.t_label = v2;
                        end
                    elseif v1 == 'u'
                        if length( v2 ) ~= self.nu
                            error( 'tsm_Data.set_Label u not string{%d}', self.nu )
                        end
                        if not( iscell( v2 ) )
                            error( 'tsm_Data.set_Label u not string{%d}', self.nu )
                        end
                        for i=1:self.nu
                            if not( ischar( v2{i} ) )
                                error( 'tsm_Data.set_Label u(%d) not string', i )
                            end
                        end
                        self.u_label = v2;
                        
                    elseif v1 == 'y'
                        if not( ischar( v2 ) )
                            error( 'tsm_Data.set_Label y not string' )
                        else
                            self.y_label = v2;
                        end
                    else
                        error( 'tsm_Data.set_Label arg 1 not t/u/y' )
                    end
                else
                    error( 'tsm_Data.set_Label not arg or arg1,arg2' )
                end
            end
        end
        
        function set_Unit( self, signal, unit )
            if not( isempty(unit) )
                if signal == 't'
                    self.t_unit = unit;
                elseif signal == 'u'
                    if self.nu == 1
                        self.u_unit{1} = unit;
                    elseif self.nu > 1 && iscell( unit )
                        for i=1:self.nu
                            self.u_unit{i} = unit{i};
                        end
                    else
                        error( 'tsm_Data/set_Unit: unit u not char(nu)')
                    end
                elseif signal == 'y'
                    self.y_unit = unit;
                end
            end
        end
        
        function set_Limit( self, signal, limit )
            if signal == 't'
                if isequal( size(limit), [1,2] )
                    self.t_limit = limit;
                else
                    error( 'tsm_Data/set_Limit: limits u not (1x2)')
                end
            elseif signal == 'u'
                if ( self.nu == 1 && isequal( size(limit), [1,2] ) ) || ...
                        ( self.nu > 1 && isequal( size(limit), [self.nu,2] ) )
                    for i=1:self.nu
                        self.u_limit(i,:) = limit(i,:);
                    end
                else
                    error( 'tsm_Data/set_Limit: limits u not (nux2)')
                end
            elseif signal == 'y'
                if isequal( size(limit), [1,2] )
                    self.y_limit = limit;
                else
                    error( 'tsm_Data/set_Limit: limits y not (1x2)')
                end
            else
                error( 'tsm_Data/set_Limit: no t/u/y limits')
            end
        end
        
    end
end