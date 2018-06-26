function crankNicolson(odefun,M,tSpan,tSteps,u0,tol,process)
%CRANKNICOLSON Solves the discrete initial value problem
%   M ( u_k - u_k-1 )/?t + 1/2*F( t_k , u_k ) + 1/2*F( t_k-1 , u_k-1 ) = 0
%   u_0 = u0
%for k = 1, 2, ..., tSteps
%   Input:
%       odefun: a function handle to a function that evaluates F(t,u) and
%       its derivative DF(t,u) with respect to u
%       M: a mass matrix
%       tSpan: a 1x2 array with the inital and final times
%       tSteps: an integer for the number of time steps such that
%                tSpan(2) - tSpan(1)
%           ?t = -------------------
%                      tSteps
%       u0: a column vector defining u_0 = u(tSpan(1))
%       tol: positive real number, tolerance for the residual at each time
%       step
%       process (optional): a function handle to three functions
%       	before(t,u)
%           during(t,u)
%           after(t,u)
%       that can be used for plotting u(t) or other commands that will be
%       executed before the first time step, at every time step and after
%       the last time step, respectively

% step size
dt = (tSpan(2)-tSpan(1))/tSteps;

if nargin > 6
    plotSnapshots = true;
    [before,during,after] = process();
    
    % preprocessing commands
    before(tSpan(1),u0);
else
    plotSnapshots = false;
end

% time stepping
for k = 1:tSteps
    % time at step k
    t = tSpan(1) + k*dt;
    
    % u at step k
    u = newton(@(u) newtonfun(t,u),u0,tol); 
    
    if plotSnapshots
        % postprocessing commands
        during(t,u);
    end
    
    % save u_k for the next iteration 
    u0 = u;
end

if plotSnapshots
    % postprocessing commands
    after(t,u);
end

    function [G,DG] = newtonfun(t,u)
        %NEWTONFUN Evalutates the function
        %   G(t,u) = M ( u - u_k-1 )/?t + F( t , u )
        %and its derivative DG(t,u) with respect to u
        [F0,DF0] = odefun(t-dt,u0);
        [F,DF] = odefun(t,u);
        G = M*(u-u0)./dt + 1/2*F+ 1/2*F0;
        DG = M./dt + 1/2*DF + 1/2*DF0;
    end
end