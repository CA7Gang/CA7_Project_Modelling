classdef PumpComponent
      %ValveComponent An object-oriented implementation of a water pump
    
    
    properties
             
        a2  % Pump flow coefficient
        a1  % Pump cross-term coefficient
        a0  % Pump velocity coefficient
        type % Component type
                               
        % Non-zero component properties for this component type (pump)
        alpha  % Pressure difference across component (only relevant for pumps)
        alpha_SI % Pressure difference across component in SI units
        
        % Obligate zero component properties for this component type (valve)
        dz  % Loss due to height change
        dz_SI  % Loss due to height change in SI units
        J % Flow dynamics term
        Lambda % Combined flow-dependent loss term
        J_SI  % J term in SI units
        Lambda_SI  % Lambda term in SI units
        mu  % Flow-dependent loss term for valves (only relevant for valves)
        mu_SI % Flow-dependent valve loss term in SI units
    end
    
    methods
        function obj = PumpComponent(a0,a1,a2)
          %PumpComponent Constructs an instance of the pump class
            if nargin > 0
            obj.a2 = a2;
            obj.a1 = a1;
            obj.a0 = a0;
            obj.type = 'Pump';
            
            % Relevant component parameters
            obj.alpha_SI = PumpModel(obj);
            obj.alpha = @(q,w) obj.alpha_SI(q,w)/(10^5*3600^2); % Inherit the function for alpha
            
            % Irrelevant component parameters
            obj.dz = 0;
            obj.dz_SI = 0;
            obj.J = 0;
            obj.J_SI = 0;
            obj.Lambda = @(q) 0;
            obj.Lambda_SI = @(q) 0;
            obj.mu = @(q,OD) 0;
            obj.mu_SI = @(q,OD) 0;
            end
            
        end
        
        function alpha = PumpModel(obj)
        %PumpModel Constructs the nonlinear model of the pump
            alpha = @(q,w) obj.a0*w^2+obj.a1*w*q+obj.a2*abs(q)*q;
        end
    end
end

