%% tsm_Conc: Definition of class for TS model premise

% $Id$

classdef tsm_Prem < tsm_Base
    
    methods
        
        %% Constructor
        function obj = tsm_Prem( model, varargin )
            % input: model

            p = inputParser;
            p.addRequired( 'model', @isa( tsModel ) )
            p.addParameter('Name','undefined',@ischar)
            p.addParameter('Comment','',@ischar)
            p.parse( model, varargin{:} )
            opts = p.Results;

            obj@tsm_Base( opts.Name, opts.Comment );            
            
            obj.setName( name );
            if nargin > 1
                obj.addComment( comment );
            end
        end
        
    end
    
end
