classdef fighandle < handle
    %FIGHANDLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fig
        params1D
        params2D
        results1D
        DA
    end
    
    methods
        
        function obj = fighandle
            obj.params1D.N   = 100;
            obj.params1D.L   = 1;
            obj.params1D.T   = 1;
            obj.params1D.y0  = 100;
            obj.params1D.xi0 = 100;

        end
    end
end

