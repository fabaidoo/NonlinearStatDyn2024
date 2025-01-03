function [d, v, ifin] = NewmarkSteps(d0, v0, T1, M, K, beta, gamma, flag)

if nargin == 7
flag = false;
end

%T1 = 6.283813688542157; %period of smallest eigenvalue of K
%T1= 0.000628255698866;
dt = T1/20; %time step for discretization%
T = 0:dt:5*T1;%1:100;%
if beta==0 && gamma==1/2
    dt = 0.01999;%0.01414;%
T = 0:100;
end

d = zeros(2, length(T));
d(:,1) = d0;
v = zeros(2, length(T));
v(:,1) = v0;
ifin = zeros(1, length(T));

d_n = d0;
v_n = v0;

for j=2:length(T)
  
    %solve for next time step
    [d_n1, v_n1, iloop] = Newmarkloop(d_n, v_n, dt,M, K, beta, gamma,...
        flag);
    
    %update history
    d(:, j) = d_n1;
    v(:, j) = v_n1;
    ifin(j) = iloop;

    %update initial condition
    d_n = d_n1;
    v_n = v_n1;
end

%{
figure
subplot(2, 1, 1)
plot(1:length(T), d(1,:), 'LineWidth',1.5)
grid
ylabel('d_1')
subplot(2, 1, 2)
plot(1:length(T), v(1,:), 'LineWidth',1.5)
grid
xlabel('n')
ylabel('v_1')

figure
subplot(2,1,1)
plot(1:length(T), d(2,:), 'LineWidth',1.5)
grid
ylabel('d_2')
subplot(2,1,2)
plot(1:length(T), v(2,:), 'LineWidth',1.5)
grid
xlabel('n')
ylabel('v_2')
%}

%ifin

end