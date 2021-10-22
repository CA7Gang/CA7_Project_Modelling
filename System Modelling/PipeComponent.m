classdef PipeComponent
    %PipeComponent An object-oriented implementation of fluid flow in a
    %pipe
    properties
        Length  % Pipe length [m]
        rho  % Fluid density [kg/m^3]
        Area  % Pipe cross-sectional area [m^2]
        Diameter  % Pipe diameter in [m]
        g  % Gravitational acceleration [m/s^2]
        f  % Flow-dependent friction factor (calculated internally)
        eta  % Pipe roughness height [m]
        Reynolds  % Flow-dependent Reynolds number of flow through pipe
        kf  % Form-loss coefficient (provided my manufacturer)
        hm  % Head loss due to form resistance
        hf  % Head loss due to surface resistance
        h  % Geodesic height change along traverse [m]
        type % Component type
        
        % Non-zero component properties for this component type (pipe)
        dz  % Loss due to height change
        dz_SI  % Loss due to height change in SI units
        J  % Flow dynamics term
        Lambda  % Combined flow-dependent loss term
        J_SI  % J term in SI units
        Lambda_SI  % Lambda term in SI units
        
        % Obligate zero component properties for this component type (pipe)
        mu % Flow-dependent loss term for valves (only relevant for valves)
        mu_SI % Flow-dependent valve loss term in SI units
        alpha % Pressure difference across component (only relevant for pumps)
        alpha_SI % Pressure difference across component in SI units 
        
        
        
    end
    
    methods
        function obj = PipeComponent(Length,rho,Diameter,g,eta,Reynolds,kf,h)
            %PipeComponent Constructs an instance of the pipe class
            
            if nargin > 0
            % Fixed object attributes
            obj.Length = Length; 
            obj.rho = rho; 
            obj.Area = pi*(Diameter/2)^2; 
            obj.Diameter = Diameter;
            obj.g = g;
            obj.eta = eta;
            obj.Reynolds = Reynolds;
            obj.kf = kf;
            obj.h = h;
            obj.type = 'Pipe';
            
            % Attributes calculated via fixed attributes
            obj.f = obj.TurbFF();
            obj.hf = obj.DarcyWeisbach(); 
            obj.hm = obj.SwameeForm();
            
            obj.Lambda_SI = obj.FlowLoss();
            obj.Lambda = @(q) obj.Lambda_SI(q)/(10^5*3600^2);
            
            obj.J_SI = obj.FlowDyn();
            obj.J = obj.J_SI/(10^5*3600^2);
            
            obj.dz_SI = obj.HeightLoss();
            obj.dz = obj.dz_SI/(10^5);
            
            % Irrelevant component properties
            obj.mu = @(q,OD) 0;
            obj.mu_SI = @(q,OD) 0;
            obj.alpha = @(q,w) 0;
            obj.alpha_SI = @(q,w) 0;
            
            end
        end
        
        function FF = TurbFF(obj)
            %TurbFF Formula for friction factor in turbulent flow regime
            FF = 1.325*(log(obj.eta/(3.7*obj.Diameter)+5.74/(obj.Reynolds^0.9)))^-2; 
        end
        
        function hf = DarcyWeisbach(obj)
            %DarcyWeissbach Calculates head loss term due to surface resistance
            %via the Darcy-Weissbach equation
            hf = obj.f*8*obj.Length/(pi^2*obj.g*obj.Diameter^5);
        end
        
        function hm = SwameeForm(obj)
            %SwameeForm Calculates head loss term due to form resistance
            %via the Swamee equation (see 2.2, e.q 2.7b in ISBN 9780470225059)
            hm = obj.kf*8/(pi^2*obj.g*obj.Diameter^4);
        end
        
        function J = FlowDyn(obj)
            %FlowDyn Calculates the flow dynamics term J
            J = obj.Length*obj.rho/obj.Area;
        end
        
        function Lambda = FlowLoss(obj)
            %FlowLoss Calculates the flow-dependent loss term Lambda
            Lambda = @(q) (obj.hf+obj.hm)*obj.rho*abs(q)*q;
        end
        
        function dz = HeightLoss(obj)
            %HeightLoss Calculates loss according to height gradient
            dz = obj.h*obj.g*obj.rho;
        end
    end
end

