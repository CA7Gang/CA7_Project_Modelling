classdef PipeComponent
    %PipeComponent An object-oriented implementation of fluid flow in a
    %pipe
    properties
        Length {mustBeNumeric} % Pipe length [m]
        rho {mustBeNumeric} % Fluid density [kg/m^3]
        Area {mustBeNumeric} % Pipe cross-sectional area [m^2]
        Diameter {mustBeNumeric} % Pipe diameter in [m]
        g {mustBeNumeric} % Gravitational acceleration [m/s^2]
        f {mustBeNumeric} % Flow-dependent friction factor (calculated internally)
        eta {mustBeNumeric} % Pipe roughness height [m]
        Reynolds {mustBeNumeric} % Flow-dependent Reynolds number of flow through pipe
        kf {mustBeNumeric} % Form-loss coefficient (provided my manufacturer)
        hm {mustBeNumeric} % Head loss due to form resistance
        hf {mustBeNumeric} % Head loss due to surface resistance
        h {mustBeNumeric} % Geodesic height change along traverse [m]
        type % Component type
        
        % Non-zero component properties for this component type (pipe)
        dz {mustBeNumeric} % Loss due to height change
        dz_SI {mustBeNumeric} % Loss due to height change in SI units
        J {mustBeNumeric} % Flow dynamics term
        Lambda {mustBeNumeric} % Combined flow-dependent loss term
        J_SI {mustBeNumeric} % J term in SI units
        Lambda_SI {mustBeNumeric} % Lambda term in SI units
        
        % Obligate zero component properties for this component type (pipe)
        mu % Flow-dependent loss term for valves (only relevant for valves)
        mu_SI % Flow-dependent valve loss term in SI units
        dp % Pressure difference across component (only relevant for pumps)
        dp_SI % Pressure difference across component in SI units 
        
        
        
    end
    
    methods
        function obj = PipeComponent(Length,rho,Area,Diameter,g,eta,Reynolds,kf,h)
            %PipeComponent Constructs an instance of the pipe class
            
            if nargin > 0
            % Fixed object attributes
            obj.Length = Length; 
            obj.rho = rho; 
            obj.Area = Area; 
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
            obj.Lambda = obj.Lambda_SI/(10^5*3600);
            
            obj.J_SI = obj.FlowDyn();
            obj.J = obj.J_SI/(10^5*3600);
            
            obj.dz_SI = obj.HeightLoss();
            obj.dz = obj.dz_SI/(10^5);
            
            % Irrelevant component properties
            obj.mu = @(OD) 0;
            obj.mu_SI = @(OD) 0;
            obj.dp = @(w) 0;
            obj.dp_SI = @(w) 0;
            
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
            Lambda = (obj.hf+obj.hm)*obj.rho;
        end
        
        function dz = HeightLoss(obj)
            %HeightLoss Calculates loss according to height gradient
            dz = obj.h*obj.g*obj.rho;
        end
    end
end

