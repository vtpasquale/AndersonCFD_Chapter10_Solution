function primitivesOut=updateBoundaryConditions(primitivesIn,inflow,Tw_Tinf)
% Update solution boundary conditions
%
% INPUTS
% primitivesIn = [Primitives] Domain primitives
% inflow       = [Primitives] Primitives at the inflow boundary
% Tw_Tinf      = [double] Wall temperature specification
%                    if Tw_Tinf > 0
%                          Tw/Tinf = Tw_Tinf
%					 else
%                          Adiabatic wall
%                    end
%
% OUTPUTS
% primitivesOut = [Primitives] Domain primitives with boundary values updated

% Anthony Ricciardi
% July 2021						   

[u,v,p,T] = primitivesIn.deal();
[Vinf,~,pinf,Tinf] = inflow.deal();

% inflow bounary, x(:,1) == 0
u(:,1) = Vinf;
v(:,1) = 0;
p(:,1) = pinf;
T(:,1) = Tinf;

% upper boundary, y(end,:) == lvert
u(end,:) = Vinf;
v(end,:) = 0;
p(end,:) = pinf;
T(end,:) = Tinf;

% outflow boundary, x(:,end) == lhori
u(:,end) = 2*u(:,end-1) - u(:,end-2);
v(:,end) = 2*v(:,end-1) - v(:,end-2);
p(:,end) = 2*p(:,end-1) - p(:,end-2);
T(:,end) = 2*T(:,end-1) - T(:,end-2);

% plate boundary, y(1,:) == 0
u(1,:) = 0;
v(1,:) = 0;
p(1,:) = 2*p(2,:) - p(3,:);

if Tw_Tinf > 0
    % constant temperature wall
    T(1,:) = Tinf * Tw_Tinf;
else 
    % adiabatic wall
    % 0 = dfdy(1,:) = (-f(3,:) + 4*f(2,:) - 3*f(1,:))./(2*dy);
    % 0 = -f(3,:) + 4*f(2,:) - 3*f(1,:) 
    % 3*f(1,:) = -f(3,:) + 4*f(2,:)
    % f(1,:) = (4*f(2,:)- f(3,:))./3
    
%     T(1,:) = (4*T(2,:)- T(3,:))./3; % [O(dy^2)]
    T(1,:) = T(2,:);  % [O(dy^2)] - Seems to match Anderson more closely
end

% leading edge, x(1,1) == y(1,1) == 0, supersedes others
u(1,1) = 0;
v(1,1) = 0;
p(1,1) = pinf;
T(1,1) = Tinf;

% Export updated primitives
primitivesOut = Primitives(u,v,p,T);