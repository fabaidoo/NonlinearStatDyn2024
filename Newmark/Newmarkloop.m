function [d_1, v_1, ifin] = Newmarkloop(d0, v0, dt, M,K,  beta,gamma, flag)

if nargin == 7
flag = false;
end

%T1 = 6.283813688542157; %period of smallest eigenvalue of K
%T= 0.000628255698866;
%dt = T1/20; %time step for discretization
if flag == false
   Mstar = M + dt^2 * beta * K;
elseif flag == true
    Mstar = M + dt^2 * beta * diag(diag(K));
end

a0 = - M \ (K * d0);

%initialize the predictor
di = d0 + dt*v0 + (dt^2/2) * (1 - 2*beta)*a0 ;
vi = v0 + dt * (1 - gamma)*a0;
ai = [0;0];

dF0 = -M*ai - K*di; %initial residual

dFi = dF0;
maxiter = 100;
tol = 1e-3;
for i = 1:maxiter
%dFi = - M*ai - K*di;

da = Mstar \ dFi; 

%corrector step
ai = ai + da;
vi = vi + gamma * dt * da;
di = di + beta * dt^2 * da;

dFi = - M*ai - K*di;

%check for convergence
if norm(dFi) < tol * norm(dF0) 
    %a_1 = ai;
    v_1 = vi;
    d_1 = di;
    ifin = i;
    break
else 
    if i == maxiter
        norm(dFi) 
        disp('Failed')
    end
    continue
end



end
















end