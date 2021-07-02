classdef Primitives
    % Class to store and interface with primitive variables for Anderson
    % Chapter 10 example. Independent primitives are selected so they are
    % consistent with those specified at the boundaries. Gas properties are
    % stored to assist with calculation of dependent properties. SI units
    % are used.
    %
    % ny and nx are the number of x and y grid points, respectively. It is
    % assumed that the 2D grid is created using the meshgrid() function.
    
	% Anthony Ricciardi
	% July 2021
	
    properties
        u % [ny,nx double] x-direction velocity
        v % [ny,nx double] y-direction velocity
        p % [ny,nx double] pressure
        T % [ny,nx double] absolute temperature
    end
    properties (Dependent = true)
        r % [ny,nx double] density
        e % [ny,nx double] internal energy
        Et % [ny,nx double] total energy
        mu % [ny,nx double] molecular viscosity coefficient
        lambda % [ny,nx double] second viscosity coefficient
        a % [ny,nx double] sound speed
        k % [ny,nx double] thermal conductivity 
    end
    properties (Constant = true)
        mu0 = 0.000017894; % [double] Pa, Sea Level Dynamic Viscosity
        T0 = 288.16; % [double] Kelvin, Sea Level Temperature
        gm = 1.4;  % [double] ratio of specific heats
        Pr = 0.71; % [double] Prandtl number
        R = 287;   % [double] J/(kg K), Specific gas constant
    end
    properties (Dependent = true)
        cv % [double] Specific Heat for an Ideal Gas at Constant volume
        cp % [double] Specific Heat for an Ideal Gas at Constant Pressure
    end
    methods % get methods for dependent properties
        function r = get.r(obj)
            % density of a perfect gas
            r = obj.p./(obj.R*obj.T);
        end
        function e = get.e(obj)
            % internal energy of a calorically perfect gas (constant specific heats)
            e = obj.cv*obj.T;
        end
        function Et = get.Et(obj)
            % total energy
            Et = obj.r.*(obj.e + .5*(obj.u.^2 + obj.v.^2));
        end
        function cv = get.cv(obj)
            % ideal gas
            cv =  obj.R/(obj.gm-1);
        end
        function cp = get.cp(obj)
            % ideal gas
            cp =  obj.gm*obj.cv;
        end
        function mu = get.mu(obj)
            % Sutherland's Law
            mu = obj.mu0.*(obj.T./obj.T0).^(3/2) .* (obj.T0 + 110)./(obj.T + 110);
        end
        function lambda = get.lambda(obj)
            % Stokes's hypothesis
            lambda = -2/3 * obj.mu;
        end
        function a = get.a(obj)
            % ideal gas sound speed
            a = sqrt(obj.gm*obj.R*obj.T); % ideal gas
        end
        function k = get.k(obj)
            % Constant Prandtl number
            k = obj.cp*obj.mu./obj.Pr;
        end
        function [muOut,lambdaOut,kOut] = getMuLambdaK(obj)
            % The dependent approach is inefficient when these three values
            % are needed at the same time. This function provides them with
            % better efficiency.
            muOut = obj.mu;
            lambdaOut = -2/3 * muOut;
            kOut = obj.cp*muOut./obj.Pr;
        end
    end
    methods
        function obj = Primitives(uIn,vIn,pIn,TIn)
            % class constructor
            obj.u = uIn; obj.v = vIn; obj.p = pIn; obj.T = TIn;
        end
        function [u,v,p,T] = deal(obj)
           % convenience function to deal variables
           u = obj.u; v = obj.v; p = obj.p; T = obj.T;
        end
        function [ReX,ReY] = calculateReynoldsNumber(obj,referenceLength)
            % ReX = [ny,nx double] Reynolds Number for x direction
            % ReY = [ny,nx double] Reynolds Number for y direction
            r_mu = obj.r./ obj.mu;
            ReX = referenceLength  .* obj.u .* r_mu;
            ReY = referenceLength  .* obj.v .* r_mu;
        end
        function dt = calculateTimeStep(obj,dx,dy,K)
            Mu = obj.mu;
            vp = max(4/3*Mu,...
                     obj.gm*Mu./obj.Pr)./obj.r; % Anderson has typos
            dtCFL = 1./( abs(obj.u)./dx + abs(obj.v)./dy + obj.a.*sqrt(1/dx^2 + 1/dy^2) + 2*vp*(1/dx^2 + 1/dy^2));
            dt = K *min(dtCFL(:));
        end
    end
end

