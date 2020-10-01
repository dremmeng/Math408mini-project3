utrue = @(t) [1+3.*exp(-8.*t); -3.*exp(-8.*t)];   % exact solution
f = @(t,u) [u(2);u(2)*(u(2)-1)/(u(1))];
U0 = [.5;-3];
tspan = [0 2];
kvals=logspace(-3,-1,10);
maxerrvec1 = [];
maxerrvec2 = [];
for i=1:length(kvals)
   k = kvals(i);
   U1 = utrue(k);
   U2 = utrue(k*2);
   % Linear multistep methods
   [tvec,Uvec1] = LMM1(f,tspan,U0,U1,k);
   [tvec,Uvec2] = LMM2(f,tspan,U0,U1,U2,k);
   utvec = utrue(tvec);   % evaluate exact solution at same points as numerical solution
   % errors:
   Evec1 = utvec - Uvec1;
   Evec2 = utvec - Uvec2;
   maxerrvec1(i) = norm(Evec1,'inf');
   maxerrvec2(i) = norm(Evec2,'inf');
   if i==1 
       fprintf('       k         LMM1 maxerror   LMM2 maxerror\n')
   end
   % print line of table:
   fprintf('%13.4e   %13.4e   %13.4e\n',k,maxerrvec1(i),maxerrvec2(i))
end
% plot errors:
loglog(kvals,maxerrvec1,kvals,maxerrvec2,'LineWidth',2)
legend('LMM1','LMM2','Location','northeast')