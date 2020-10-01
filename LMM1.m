function [tvec,Uvec] = LMM1(f,tspan,U0,U1,k)
% [tvec,Uvec] = LMM1(f,tspan,U0,U1,k)
% Inputs
% f: name or function handle of the right-hand side function f:(t,u)->f(t,u)
% tspan(1),U0: initial condition (U0 can be s-vector)
% U1: additional starting value from exact solution
% tspan(2): end time, so that number of steps N = (tspan(2)-tspan(1))/k
% k: stepsize
% Outputs
% tvec: vector of t values
% Uvec: vector (or matrix) of corresponding u values

U0 = U0(:);          % make sure U0 is a column vector
s = length(U0);      % number of equations in system
tvec = tspan(1):k:tspan(2);   % a row vector
N = length(tvec)-1;
Uvec = zeros(s,N+1);
% start up with exact values (provided as input)
Uvec(:,1) = U0;
Uvec(:,2) = U1;
for n = 1:N-1
    fvalue = f(tvec(n),Uvec(:,n))
    fvalu = f(tvec(n+1),Uvec(:,n+1))
    Uvec(:,n+2) = Uvec(:,n+1)+ (k/3)*(-2*fvalue + 3*fvalu)
end
tvec = tvec';
Uvec = Uvec';        % to match MATLAB output