function NonlinearMethodPlotter(x)
%returns various plots for the problem

Fi = linspace(0, 10, 41);

d_Exact = zeros(1, length(Fi));
d_NR_ls = zeros(1, length(Fi));
d_NR_nls = zeros(1, length(Fi));
d_MNR_ls = zeros(1, length(Fi));
d_MNR_nls = zeros(1, length(Fi));
d_BFGS_ls = zeros(1, length(Fi));
d_BFGS_nls = zeros(1, length(Fi));

i_NR_ls = zeros(1, length(Fi));
i_NR_nls = zeros(1, length(Fi));
i_MNR_ls = zeros(1, length(Fi));
i_MNR_nls = zeros(1, length(Fi));
i_BFGS_ls = zeros(1, length(Fi));
i_BFGS_nls = zeros(1, length(Fi));

d_Exact(1) = ExactN(x, Fi(1));
[d_NR_ls(1), i_NR_ls(1), ~] = NewtonRaphLineSearch(Fi(1), 1, x );
[d_NR_nls(1), i_NR_nls(1), ~] = NewtonRaphNoLineSearch(Fi(1), 1, x );
[d_MNR_ls(1), i_MNR_ls(1), ~] = ModifiedNewtonRaphLineSearch(Fi(1), 1, x);
[d_MNR_nls(1), i_MNR_nls(1), ~] = ModifiedNewtonRaphNoLineSearch(Fi(1), 1, x);
[d_BFGS_ls(1), i_BFGS_ls(1), ~] = MNRwithBFGSlineSearchable(Fi(1), 1, x, true);
[d_BFGS_nls(1), i_BFGS_nls(1), ~] = MNRwithBFGSlineSearchable(Fi(1), 1, x, false);


for k = 2: length(Fi)
    d_Exact(k) = ExactN(x, Fi(k));
    [d_NR_ls(k), i_NR_ls(k), ~] ...
        = NewtonRaphLineSearch(Fi(k),d_NR_ls(k-1), x );

    [d_NR_nls(k), i_NR_nls(k), ~] ...
        = NewtonRaphNoLineSearch(Fi(k), d_NR_nls(k-1), x );

    [d_MNR_ls(k), i_MNR_ls(k), ~] ...
        = ModifiedNewtonRaphLineSearch(Fi(k), d_MNR_ls(k-1), x);

    [d_MNR_nls(k), i_MNR_nls(k), ~] ...
        = ModifiedNewtonRaphNoLineSearch(Fi(k), d_MNR_nls(k-1), x);

    [d_BFGS_ls(k), i_BFGS_ls(k), ~] ...
        = MNRwithBFGSlineSearchable(Fi(k), d_BFGS_ls(k-1), x, true);

    [d_BFGS_nls(k), i_BFGS_nls(k), ~] ...
        = MNRwithBFGSlineSearchable(Fi(k), d_BFGS_ls(k-1), x, false);

end

NX = NFunc(d_Exact);
N_NR_ls = NFunc(d_NR_ls);
N_NR_nls = NFunc(d_NR_nls);
N_MNR_ls = NFunc(d_MNR_ls);
N_MNR_nls = NFunc(d_NR_ls);
N_BFGS_ls = NFunc( d_BFGS_ls);
N_BFGS_nls = NFunc(d_BFGS_nls);







figure(Name='Calculated Loads')
subplot(1, 2, 1)
plot(d_Exact, NX,  'k.', 'MarkerSize', 15, 'DisplayName', 'Exact')
hold on

plot(d_NR_ls, N_NR_ls, '.-', 'DisplayName', 'N-R LS', 'LineWidth',1.5,...
    'Color',[0 0.4470 0.7410])
%
plot(d_NR_nls, N_NR_nls, 'o--', 'DisplayName', 'N-R no LS', 'LineWidth',1.5,...
    'Color',[0.8500 0.3250 0.0980])

plot(d_MNR_ls, N_MNR_ls, '*-', 'DisplayName', 'MN-R LS', 'LineWidth',1.5,...
    'Color',[0.9290 0.6940 0.1250])
plot(d_MNR_nls, N_MNR_nls, 'x--', 'DisplayName', 'MN-R no LS', 'LineWidth',1.5,...
    'Color',[0.4940 0.1840 0.5560])

plot(d_BFGS_ls, N_BFGS_ls, 'd-', 'DisplayName', 'BFGS LS', 'LineWidth',1.5,...
    'Color',[0.4660 0.6740 0.1880])
plot(d_BFGS_nls, N_BFGS_nls, '^--', 'DisplayName', 'BFGS no LS', 'LineWidth',1.5,...
    'Color',[0.3010 0.7450 0.9330])
%}
title('Calculated d')
grid on
xlabel('d')
ylabel('N(d)')
legend('FontSize', 15)

subplot(1,2,2)

semilogy(d_NR_ls, abs(N_NR_ls - NX)./abs(NX), '.', 'DisplayName', 'N-R LS','MarkerSize', 10,...
    'Color',[0 0.4470 0.7410], 'LineWidth', 1.5)
hold on
semilogy(d_NR_nls, abs(N_NR_nls - NX)./abs(NX), 'o', 'DisplayName', 'N-R no LS','MarkerSize', 10,...
    'Color',[0.8500 0.3250 0.0980], 'LineWidth', 1.5)

semilogy(d_MNR_ls,abs(N_MNR_ls-NX)./abs(NX), '*', 'DisplayName', 'MN-R LS','MarkerSize', 10,...
    'Color',[0.9290 0.6940 0.1250], 'LineWidth', 1.5)
semilogy(d_MNR_nls, abs(N_MNR_nls-NX)./abs(NX), 'x', 'DisplayName', 'MN-R no LS','MarkerSize', 10,...
    'Color',[0.4940 0.1840 0.5560], 'LineWidth', 1.5)

semilogy(d_BFGS_ls, abs(N_BFGS_ls -NX)./abs(NX), 'd', 'DisplayName', 'BFGS LS','MarkerSize', 10,...
    'Color',[0.4660 0.6740 0.1880], 'LineWidth', 1.5)
semilogy(d_BFGS_nls, abs(N_BFGS_nls-NX)./abs(NX), '^', 'DisplayName', 'BFGS no LS','MarkerSize', 10,...
    'Color',[0.3010 0.7450 0.9330], 'LineWidth', 1.5)
grid on
title(sprintf('Relative Error x = %i', x))
xlabel('d')
ylabel('Relative Error')

sgtitle(sprintf('x = %i', x))


figure(Name='Convergence')


plot(Fi, i_NR_ls, '.-', 'DisplayName', 'N-R LS', 'LineWidth',1.5,...
    'Color',[0 0.4470 0.7410])
hold on
plot(Fi, i_NR_nls, 'o--', 'DisplayName', 'N-R no LS', 'LineWidth',1.5,...
    'Color',[0.8500 0.3250 0.0980])

plot(Fi, i_MNR_ls, '*-', 'DisplayName', 'MN-R LS', 'LineWidth',1.5,...
    'Color',[0.9290 0.6940 0.1250])
plot(Fi, i_MNR_nls, 'x--', 'DisplayName', 'MN-R no LS', 'LineWidth',1.5,...
    'Color',[0.4940 0.1840 0.5560])

plot(Fi, i_BFGS_ls, 'd-', 'DisplayName', 'BFGS LS', 'LineWidth',1.5,...
    'Color',[0.4660 0.6740 0.1880])
plot(Fi, i_BFGS_nls, '^--', 'DisplayName', 'BFGS no LS', 'LineWidth',1.5,...
    'Color',[0.3010 0.7450 0.9330])
title(sprintf('Number of iterations vs load step, x = %i', x))
legend('FontSize', 15)
grid on


    function Nd = NFunc(di)
        p = [0.5, -5, (x+Fi), -10*Fi];
        Nd  = polyval(p, di);

       

    end











end