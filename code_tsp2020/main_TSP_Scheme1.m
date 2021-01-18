%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% code: CS-based Active User Detection and Channel Estimation for Massive Access (Scheme 1)
% spatial-frequency structured sparsity
% written by Malong Ke (kemalong@bit.edu.cn), Beijing Institute of Technology
% version: 2019.03.14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Prepare the workspace && Set parameters 
clc; clear; close all
% parameters of system && channel model
K = 500;           % number of potential UDs in the cell
Ka = 50;           % number of active UDs
N_BS = 16;         % number of BS antennas
P = 1;             % number of pilot subcarriers
SNR = 30;          % SNR in dB
T_max = K;         % maximum time slot overhead
N = 2048;          % number of subcarriers
Bw = 10e6;         % system bandwidth
fs = Bw;           % frequency of sampling
L_cp = 64;         % length of cyclic prefix
sigma2_alpha = 1;  % average power of path, generally, sigma2_alpha = 1;
Lp_min = 8;        % number of paths
Lp_max = 14;
% parameters of GMMV-AMP algorithm
damp = 0.3;       % control the convergence speed (prevent algorithm from diverging)
sel_NN_set = 0;   % select nearest neighbor set for NNSPL 0: strctured sparsity   1: clustered sparsity
niter = 200;      % maximum number of AMP iterations
tol = 1e-5;       % termination threshold
% parameters of AUD
p_cg = 0.9;
p_bi = 0.9;
epsilon_bi = 0.5;
% parameters of SE
K_mc = 1000;      % number of Monte Carlo samples for SE

%% Simulations
tic
N_sim = 10;
T_set = 30:5:80;                   % set time slot overhead
Pe_bi = zeros(numel(T_set),1);     % EDP of BI-AD
Pe_cg = zeros(numel(T_set),1);     % EDP of CG-AD
MSE_amp = zeros(numel(T_set),1);  
for sim_ii = 1:N_sim
   %% Generate simulation samples
    % generate channel model
    act_flag = zeros(K,1);       
    index = randperm(K);
    act_flag(index(1:Ka)) = 1; 
    H = zeros(K,N_BS,P);
    for k = 1:K 
        Lp = randi([Lp_min,Lp_max],1,1);
        [H_f, H_ang_f] = generate_FSF_channel_grid_match(N_BS, N, P, Bw, fs, L_cp, sigma2_alpha, Lp);
        H(k,:,:) = act_flag(k).*H_f;
    end
    % pilot matrix 
    gamma = 1;
    S = sqrt(gamma/2).*(randn(T_max,K,P)+1i.*randn(T_max,K,P));
    % noisy received signals
    Y = zeros(T_max,N_BS,P);
    for p = 1:P
        Y(:,:,p) = awgn(S(:,:,p)*H(:,:,p), SNR, 'measured');
    end
     
   %% Active User Detection and Channel Estimation
    for T_ii = 1:numel(T_set)
        T = T_set(T_ii);
        
        % GMMV-AMP
        [H_amp, lambda] = gmmv_amp(Y(1:T,:,:), S(1:T,:,:), damp, niter, tol, sel_NN_set);
        % BI-AD
        act_bi = zeros(K,1);
        act_bi(sum(sum(lambda>epsilon_bi,3),2)./N_BS./P >= p_bi) = 1;
        % CG-AD
        act_cg = zeros(K,1);
        epsilon_cg = 0.01*max(max(max(abs(H_amp))));
        act_cg(sum(sum(abs(H_amp)>epsilon_cg,3),2)./N_BS./P >= p_cg) = 1;        
                                
        
        % Data recording
        Pe_bi(T_ii) = Pe_bi(T_ii) + sum(abs(act_bi-act_flag))/K;
        Pe_cg(T_ii) = Pe_cg(T_ii) + sum(abs(act_cg-act_flag))/K;
        MSE_amp(T_ii) = MSE_amp(T_ii) + norm(H_amp(1:end)-H(1:end))^2/(K*N_BS*P);
        
        fprintf('sim = %d, T = %d, Pe_bi = %4.5f, Pe_cg = %4.5f, MSE_amp = %7.9f\n', sim_ii, T, ...
                sum(abs(act_bi-act_flag))/K, sum(abs(act_cg-act_flag))/K, norm(H_amp(1:end)-H(1:end))^2/(K*N_BS*P));
        
    end 
end
toc

save Data_Scheme1_M64_P1.mat K Ka N_BS P SNR damp epsilon_bi K_mc Lp_min ...
     Lp_max sel_NN_set T_set Pe_bi Pe_cg MSE_amp N_sim

Pe_bi = Pe_bi./N_sim;
Pe_cg = Pe_cg./N_sim;
MSE_amp = MSE_amp./N_sim;

% figure Pe
figure
semilogy(T_set, Pe_bi, 'r--x', 'linewidth', 1.5);
hold on
semilogy(T_set, Pe_cg, 'b-x', 'linewidth', 1.5);
set(gca, 'GridLineStyle', '-.');
grid on;
xlabel('Time Overhead');ylabel('P_e');
legend('BI-AD', 'CG-AD');

% figure MSE
figure
semilogy(T_set, MSE_amp, 'b--o', 'linewidth', 1.5);
set(gca, 'GridLineStyle', '-.');
grid on;
xlabel('Time Overhead');ylabel('MSE');
legend('AMP');





















