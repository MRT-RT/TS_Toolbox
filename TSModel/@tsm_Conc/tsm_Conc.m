%% tsm_Conc: Definition of class for TS model conclusion

% $Id$

classdef tsm_Conc < tsm_Base
   
    properties
        model        % TS model 
    end
    
    methods
        
        %% Constructor
        function obj = tsm_Conc( model, varargin )
            % input: model

            isTSModel = @(x) isa(x,'TSModel');
            p = inputParser;
            p.addRequired( 'model', isTSModel )
            p.addParameter('Name','premise part',@ischar)
            p.addParameter('Comment','',@ischar)
            p.parse( model, varargin{:} )
            opts = p.Results;

            obj@tsm_Base( opts.Name, opts.Comment );            
            obj.model = opts.model;
        end
        
    end
    
end
