classdef ValveComponent
    %ValveComponent An object-oriented implementation of fluid flow through
    %a valve
    
    properties
        
        type % Component type
        valve_type % String denoting valve model type. Valid values are "EP" and "Linear"
        params {mustBeNumeric} % Vector containing the relevant parameters for the assigned valve type.
                               % Please specify as [tau,Vmax] if valve is EP
                               
        % Non-zero component properties for this component type (valve)
        mu  % Flow-dependent loss term for valves (only relevant for valves)
        mu_SI % Flow-dependent valve loss term in SI units
        
        % Obligate zero component properties for this component type (valve)
        dz {mustBeNumeric} % Loss due to height change
        dz_SI {mustBeNumeric} % Loss due to height change in SI units
        J {mustBeNumeric} % Flow dynamics term
        Lambda {mustBeNumeric} % Combined flow-dependent loss term
        J_SI {mustBeNumeric} % J term in SI units
        Lambda_SI {mustBeNumeric} % Lambda term in SI units
        dp % Pressure difference across component (only relevant for pumps)
        dp_SI % Pressure difference across component in SI units 
    end
    
    methods
        function obj = ValveComponent(type,params)
            %ValveComponent Constructs an instance of the valve class
            if nargin > 0
            obj.valve_type = type;
            obj.params = params;
            obj.type = 'Valve';
            
            % Relevant component parameters
            obj.mu_SI = ValveTypeSelector(obj);
            obj.mu = @(OD) obj.mu_SI(OD)/(10^5*3600); % Inherit whichever function was selected for mu
            
            % Irrelevant component parameters
            obj.dz = 0;
            obj.dz_SI = 0;
            obj.J = 0;
            obj.J_SI = 0;
            obj.Lambda = 0;
            obj.Lambda_SI = 0;
            obj.dp = @(w) 0;
            obj.dp_SI = @(w) 0;
            
               
                
            end
        end
        
        function mu = ValveTypeSelector(obj)
            %ValveTypeSelector Specifies the valve model type by keyword
             if strcmp(obj.valve_type,'EP') == 1 % Equal-percentage valve type
                 if size(obj.params) < 2
                     error('Too few parameters specified for EP valve')
                 end
                    mu = @(OD) exp(log(obj.params(1)*OD))/obj.params(1)*obj.params(2);
                elseif strcmp(obj.valve_type,'Linear') == 1 % Linear valve type
                    if size(obj.params,1) > 1
                        warning('Too many parameters specified for linear valve')
                    end
                    mu = @(OD) (1/(obj.params(1)*OD)^2);
             else
                    warning('No valve type specified. Defaulting to linear valve')    
                    mu = @(OD) (1/(obj.params(1)*OD)^2);
             end
        end
    end
end

