classdef PumpComponent
      %ValveComponent An object-oriented implementation of a water pump
    
    
    properties
             
        a2 {mustBeNumeric} % Pump flow coefficient
        a0 {mustBeNumeric} % Pump velocity coefficient
                               
        % Non-zero component properties for this component type (pump)
        dp  % Pressure difference across component (only relevant for pumps)
        dp_SI % Pressure difference across component in SI units
        
        % Obligate zero component properties for this component type (valve)
        dz {mustBeNumeric} % Loss due to height change
        dz_SI {mustBeNumeric} % Loss due to height change in SI units
        J {mustBeNumeric} % Flow dynamics term
        Lambda {mustBeNumeric} % Combined flow-dependent loss term
        J_SI {mustBeNumeric} % J term in SI units
        Lambda_SI {mustBeNumeric} % Lambda term in SI units
        mu  % Flow-dependent loss term for valves (only relevant for valves)
        mu_SI % Flow-dependent valve loss term in SI units
    end
    
    methods
        function obj = PumpComponent(a2,a0)
          %PumpComponent Constructs an instance of the pump class
            if nargin > 0
            obj.a2 = a2;
            obj.a0 = a0;
            
            % Relevant component parameters
            obj.dp_SI = PumpModel(obj);
            obj.dp = @(w) obj.dp_SI(w)/(10^5*3600); % Inherit the function for dp
            
            % Irrelevant component parameters
            obj.dz = 0;
            obj.dz_SI = 0;
            obj.J = 0;
            obj.J_SI = 0;
            obj.Lambda = 0;
            obj.Lambda_SI = 0;
            obj.mu = @(OD) 0;
            obj.mu_SI = @(OD) 0;
            end
            
        end
        
        function dp = PumpModel(obj)
        %PumpModel Constructs the nonlinear model of the pump
            dp = @(w) obj.a0*w^2;
        end
    end
end

