function [tvec,Uvec] = RK2Sys(f,tspan,U0,k)
% [tvec,Uvec] = RK2Sys(f,tspan,U0,k)
% Second-order midpoint Runge-Kutta method
% Inputs
% f: name or function handle of the right-hand side function f:(t,u)->f(t,u)
% tspan(1),U0: initial condition (u0 can be s-vector)
% tspan(2): end time, so that number of steps N = (tspan(2)-tspan(1))/k
% k: stepsize
% Outputs
% tvec: vector of t values
% Uvec: vector (or matrix) of corresponding U values

U0 = U0(:);          % make sure U0 is a column vector
s = length(U0);      % number of equations in system
tvec = tspan(1):k:tspan(2);   % a row vector
N = length(tvec)-1;
Uvec = zeros(s,N+1);
Uvec(:,1) = U0;
for n=1:N
   V1 = f(tvec(n), Uvec(:,n)); V1 = V1(:);   % force into column vectors
   V2 = f(tvec(n)+k/2, Uvec(:,n) + k/2*V1); V2 = V2(:);
   Uvec(:,n+1) = Uvec(:,n)+k*V2;
end
tvec = tvec';
Uvec = Uvec';        % to match MATLAB output

