%% tsm_Conc: Definition of class for TS model premise

% $Id$

classdef tsm_Prem < tsm_Base

    properties
        model        % TS model 
    end
    
    methods
        
        %% Constructor
        function obj = tsm_Prem( model, varargin )
            % input: model
            
            isTSModel = @(x) isa(x,'TSModel');
            p = inputParser;
            p.addRequired( 'model', isTSModel )
            p.addParameter('Name','undefined',@ischar)
            p.addParameter('Comment','',@ischar)
            p.parse( model, varargin{:} )
            opts = p.Results;
            
            obj@tsm_Base( opts.Name, opts.Comment );
            obj.model = opts.model;
            
        end
        
    end
    
end
