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
        dz {mustBeNumeric} % Geodesic height change along pipe traverse [m]
        J {mustBeNumeric} % Flow dynamics term
        
        
    end
    
    methods
        function obj = PipeComponent(Length,rho,Area,Diameter,g,eta,Reynolds,kf,dz)
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
            obj.dz = dz;
            
            % Attributes calculated via fixed attributes
            obj.f = obj.TurbFF();
            obj.hf = obj.DarcyWeisbach(); 
            obj.hm = obj.SwameeForm();
            obj.J = obj.FlowDyn();
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
        
        function dQ = calcdQ(obj,q,dP)
            %CalcdQSI Returns the change in flow with respect to time in SI
            %units
            lambda = ((obj.hf+obj.hm)*obj.rho)/(10e5*3600);
            zeta = (obj.dz*obj.g*obj.rho)/(10e5);
            J = obj.J/(10e5*3600);
            dQ = (dP-(lambda*abs(q)*q)-zeta)/(J);
        end
        
        function dQSI = calcdQSI(obj,q,dP)
            %CalcdQSI Returns the change in flow with respect to time in SI
            %units
            lambda = (obj.hf+obj.hm)*obj.rho;
            zeta = obj.dz*obj.g*obj.rho;
            dQSI = (dP-(lambda*abs(q)*q)-zeta)/obj.J;
        end
    end
end

