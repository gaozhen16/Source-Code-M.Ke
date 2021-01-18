clc; clear; close all

% detection error probability
figure 
semilogy(T_set,Perform.Pe,'k-*', 'linewidth', 1.2);
hold on
semilogy(T_set,Perform.Pe_sic,'k-.o', 'linewidth', 1.2);
hold on
set(gca, 'GridLineStyle', '-.');
grid on;
xlabel('Time Overhead $G$'); ylabel('\itP_e')

semilogy(T_set,Perform.Pe,'r-*', 'linewidth', 1.2);
hold on
semilogy(T_set,Perform.Pe_sic,'r-.o', 'linewidth', 1.2);


% NMSE
figure 
plot(T_set,Perform.NMSE,'k-o',  T_set,Perform.NMSE_ang,'r-o', T_set,Perform.NMSE_ang_sic,'g-o', 'linewidth', 1.2);
hold on
set(gca, 'GridLineStyle', '-.');
grid on;
xlabel('Time Overhead $G$'); ylabel('NMSE')

plot(T_set,Perform.NMSE,'k-*',  T_set,Perform.NMSE_ang,'r-*', T_set,Perform.NMSE_ang_sic,'g-*', 'linewidth', 1.2);


%% ref
% detection error probability
figure 
semilogy(T_set,Perform.Pe,'k-o', T_set,Perform.PM,'k--o', T_set,Perform.FA,'k-.o', 'linewidth', 1.2);
hold on
semilogy(T_set,Perform.Pe_sic,'r-x', T_set,Perform.PM_sic,'r--x', T_set,Perform.FA_sic,'r-.x', 'linewidth', 1.2);
hold on
set(gca, 'GridLineStyle', '-.');
grid on;
xlabel('Time Overhead $G$'); ylabel('\itP_e')

% NMSE
figure 
plot(T_set,Perform.NMSE,'k-o',  T_set,Perform.NMSE_ang,'r-o', T_set,Perform.NMSE_ang_sic,'g-o', 'linewidth', 1.2);
hold on
set(gca, 'GridLineStyle', '-.');
grid on;
xlabel('Time Overhead $G$'); ylabel('NMSE')
