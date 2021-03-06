classdef ValveComponent
    %ValveComponent An object-oriented implementation of fluid flow through
    %a valve
    
    properties
        
        type % Component type
        valve_type % String denoting valve model type. Valid values are "EP" and "Linear"
        params  % Vector containing the relevant parameters for the assigned valve type.
                               % Please specify as [tau,Vmax] if valve is EP
                               
        % Non-zero component properties for this component type (valve)
        mu  % Flow-dependent loss term for valves (only relevant for valves)
        mu_SI % Flow-dependent valve loss term in SI units
        
        % Obligate zero component properties for this component type (valve)
        dz  % Loss due to height change
        dz_SI  % Loss due to height change in SI units
        J  % Flow dynamics term
        Lambda  % Combined flow-dependent loss term
        J_SI  % J term in SI units
        Lambda_SI  % Lambda term in SI units
        alpha % Pressure difference across component (only relevant for pumps)
        alpha_SI % Pressure difference across component in SI units 
    end
    
    methods
        function obj = ValveComponent(type,params)
            %ValveComponent Constructs an instance of the valve class
            if nargin > 0
            obj.valve_type = type;
            obj.params = params;
            obj.type = 'Valve';
            
            % Relevant component parameters
            obj.mu = ValveTypeSelector(obj);
            obj.mu_SI = @(q,OD) obj.mu(q,OD)/(10^5*3600^2); % Inherit whichever function was selected for mu
            
            % Irrelevant component parameters
            obj.dz = 0;
            obj.dz_SI = 0;
            obj.J = 0;
            obj.J_SI = 0;
            obj.Lambda = @(q) 0;
            obj.Lambda_SI = @(q) 0;
            obj.alpha = @(q,w) 0;
            obj.alpha_SI = @(q,w) 0;
            
               
                
            end
        end
        
        function mu = ValveTypeSelector(obj)
            %ValveTypeSelector Specifies the valve model type by keyword
             if strcmp(obj.valve_type,'EP') == 1 % Equal-percentage valve type
                 if size(obj.params) < 2
                     error('Too few parameters specified for EP valve')
                 end
                    mu = @(q,OD) exp(log(obj.params(1)*OD))/obj.params(1)*obj.params(2)*abs(q)*q;
                elseif strcmp(obj.valve_type,'Linear') == 1 % Linear valve type
                    if size(obj.params,1) > 1
                        warning('Too many parameters specified for linear valve')
                    end
                    mu = @(q,OD) (1/(obj.params(1)*OD)^2)*abs(q)*q;
             else
                    warning('No valve type specified. Defaulting to linear valve')    
                    mu = @(q,OD) (1/(obj.params(1)*OD)^2)*abs(q)*q;
             end
        end
    end
end

