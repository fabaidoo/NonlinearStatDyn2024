function NewmarkImplement17

m1 = 1;
m2 = 1;
M = [m1 0; 0 m2];%mass matrix

%{
k2 = 10^4;
k1=1;
K = [k1+k2 -k2; -k2 k2];
%}
k1 = 10^4;
k2=1;
K = [k1+k2 -k2; -k2 k2];

omega = eig(K);
omega1 = omega(1);
T1 = 2*pi/sqrt(omega1); %period of smallest eigenvalue of K
dt = T1/20; %time step for discretization%
T = 0:dt:5*T1;

d0 = [10;11];
v0 = [0;0];

Method = ["Central Differences", "Trapezoidal Rule", "Damped Newmark"];
betas = [0, 1/4, 0.3025];
gammas = [1/2, 1/2, 0.6];

for i = 1: length(Method)
    [d, v, ifin] = NewmarkSteps(d0, v0, T1, M, K, betas(i), gammas(i));

    figure('Name', Method(i))
    sgtitle(Method(i))
    subplot(2, 2, 1)
    plot(1:length(T), d(1,:), 'LineWidth',1.5)
    grid
    ylabel('d_1')
    subplot(2, 2, 2)
    plot(1:length(T), v(1,:), 'LineWidth',1.5)
    grid
    xlabel('n', 'FontSize', 15)
    ylabel('v_1')

    %figure('Name',Method(i))
    %sgtitle(Method(i))
    subplot(2,2,3)
    plot(1:length(T), d(2,:), 'LineWidth',1.5)
    grid
    ylabel('d_2')
    subplot(2,2,4)
    plot(1:length(T), v(2,:), 'LineWidth',1.5)
    grid
    xlabel('n', 'FontSize', 15)
    ylabel('v_2')

    figure('Name', Method(i))
    plot(0:(length(T)-1), ifin, 'LineWidth', 1.5)
    ylabel('Number of iterations')
    xlabel('n', 'FontSize', 15)
    grid
end


end