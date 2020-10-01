function [tvec,Uvec] = LMM2(f,tspan,U0,U1,U2,k)
% [tvec,Uvec] = LMM2(f,tspan,U0,U1,U2,k)
% Inputs
% f: name or function handle of the right-hand side function f:(t,u)->f(t,u)
% tspan(1),U0: initial condition (U0 can be s-vector)
% U1,U2: additional starting values from exact solution
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
Uvec(:,3) = U2;
W = Uvec(:,3);
tol = 1e-3;
for n= 2:N
    converged = 1;
    iter = 1;
    while (converged > tol) && (iter < 10)
        fval = f(tvec(n-1),Uvec(:,n-1));
        fvalu = f(tvec(n),Uvec(:,n));
        fvalue = f(tvec(n+1), W);
        W = -1*Uvec(:,n) + 2*Uvec(:,n-1) + (k/4)*(3*fval+8*fvalu+fvalue);
        if iter > 1
            converged = norm(Wold - W, 'inf');
        end
        Wold = W;
        iter = iter +1;
    end
    Uvec(:,n+1) = W;
end
tvec = tvec';
Uvec = Uvec';        % to match MATLAB output