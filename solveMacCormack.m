function [primitives,massDiffCheck,converged,i,rumTime] = ...
          solveMacCormack(primitives,inflow,Tw_Tinf,K,x,y,maxiter)
% Solve Navier-Stokes equations for supersonic flow over a flat plate using MacCormack method.
%
% INPUTS
% primitives = [Primitives] Initial domain primitives
% inflow     = [Primitives] Primitives at the inflow boundary
% Tw_Tinf    = [double] Wall temperature specification
%                  if Tw_Tinf > 0
%                        Tw/Tinf = Tw_Tinf
%				 else
%                        Adiabatic wall
%                  end
% K          = [double] Courant number (0.5 <= K <= 0.8 recommended)
% x          = [ny,nx double] Grid point x locations (must be uniform spacing)
% x          = [ny,nx double] Grid point y locations (must be uniform spacing)
% maxiter    = [double] maximum number of time steps
% 
% OUTPUTS
% primitives    = [Primitives] Final domain primitives
% massDiffCheck = [double] Percent difference between inflow and outflow mass
% converged     = [logical] True if convergence criteria satisfied
% i             = [double] number of iterations
% rumTime       = [double] solution runtime in seconds
%
% ny and nx are the number of x and y grid points, respectively. It is
% assumed that the 2D grid is created using the meshgrid() function.

% Anthony Ricciardi
% July 2021						   
tic

% mesh data
[ny,nx] = size(x);
dx = x(1,2)-x(1,1);
dy = y(2,1)-y(1,1);
inX = 2:nx-1; % interior x
inY = 2:ny-1; % interior y

% set boundary conditions
primitives=updateBoundaryConditions(primitives,inflow,Tw_Tinf);

% time march
i = 0;
converged = false;
while ~converged && i < maxiter
    i = i + 1;
    
    % time step size
    dt = primitives.calculateTimeStep(dx,dy,K);
    
    % solution vector
    U = calculateU(primitives);
    
    % flux vectors
    E = calculateE(primitives,dx,dy,'backward');
    F = calculateF(primitives,dx,dy,'backward');
    
    % forward finite difference predictor
    dUdt_predictor = -ddxf(E,dx) - ddyf(F,dy);
    
    % Predictor step and corrector vector calculations
    U2 = U;
    U2(inY,inX,:) = U(inY,inX,:) + dt*dUdt_predictor(inY,inX,:); % interior points only
    primitives2=decodeSolutionVector(U2);
    E2 = calculateE(primitives2,dx,dy,'forward');
    F2 = calculateF(primitives2,dx,dy,'forward');
    
    % backward finite difference corrector
    dUdt_corrector = -ddxb(E2,dx) - ddyb(F2,dy);
    
    % MacCormack solution step
    dUdt = 0.5*(dUdt_predictor + dUdt_corrector);
    U(inY,inX,:) = U(inY,inX,:) + dt*dUdt(inY,inX,:);% interior points only
    primitives=decodeSolutionVector(U);
    primitives=updateBoundaryConditions(primitives,inflow,Tw_Tinf);
    
    % check density convergence
    rCurrent = primitives.r;
    if i > 1
        deltaR = max(max(abs(rCurrent - rLast)));
        if deltaR < 1e-8
            converged = true;
        end
        fprintf(1,'Iteration:%5d | delta rho: %8e\n',i,deltaR);
    end
    rLast = rCurrent;
    
end
runTime = toc;

% Mass Flow Check
massIn = trapz(y(:,1),primitives.u(:,1).*primitives.r(:,1));
massOut = trapz(y(:,end),primitives.u(:,end).*primitives.r(:,end));
massDiffCheck = 100*abs(massIn-massOut)./massIn;
fprintf(1,'Mass inflow matches mass outflow within %.3f%%.\n',massDiffCheck);
fprintf(1,'Runtime: %.2f seconds.\n',runTime);