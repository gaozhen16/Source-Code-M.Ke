%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% active user etection and channel estimation for massive access scenarios
% cell-free massive MIMO, considering large scale fading;
% exploit the sporadic traffic of users and the virtual angular domian
% sparsity of massive MIMO channels (consider power leakage caused by off-grid);
% OFDM is employed to combat time dispersive channels

% version: 2020.01.10
% by Malong Ke (kemalong@bit.edu.cn), Beijing Institute of Technology

% run time:  Nsim - 100; Worker - ; Time - 10h; Initial Time - 20/01/03 00:00 
% Nsim = 200; Kc = 400; N_BS = 16; P = 1; Time = 5h;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all
tic

%% System Parameters
Kc = 400;        % number of potential user equipments (UEs) in each small cell
Sa = 0.05;       % sparsity ratio of UE activity
N_BS = 16;       % number of access point (AP) antennas
N_UE = 1;        % number of UE antennas
% network and cell 
layer = 2;       % number of small cell layers
r_max = 1;       % radius of small cell in km
r_min = 0.01;    % minimum user-AP distance 
% path loss
loss_a = 128.1;  % 3GPP standard
loss_b = 37.6;
% transmit power and background noise
p_ad = 1;        % adaptive power p_ad = 0, 1, 2
PA = 1e3;        % power amplification
Ps_dBm = 23;     % maximum transmit power in dBm
Psd_n = -174;    % AWGN Psd dBm/Hz
% channel model
P = 1;           % nmeber of pilot subcarriers (P = N_CP)
N = 2048;        % number of subcarriers
N_CP = 64;       % cyclic prefix length
Bw = 10e6;       % bandwidth     
Lp_min = 40;     % number of MPCs
Lp_max = 100;
ang_spread = 10; % angular spread in degree
var_gain = 1;    % variance of path gain CN(beta;0,var_gain)
rho_rdd_r = 1;   % parameter of redundant distionary   
rho_rdd_t = 1;
% AMP algorithm and activity detector
Nco = 3;         % number of APs for cooperation (disable)
niter_amp = 200; % maximum number of iterations for AMP
damp_aud = 0.3;  % damping
damp_ce = 0.3;
nns_sel_aud = 1; % nearest neighbor set for refining sparsity ratio 0-each AP opearte seperately 1-Cloud  2-Edge_new  3-Edge  4-clustered
nns_sel_ce = 4;
tol_iter = 1e-5; % stop criteria
tol_re = 1e-4;
th_bi = 0.1;     % threshold of activity detector
th_rel = 0.9;    % threshold to decide reliable active user set
% SIC 
niter_sic = 3;   % maximum number of iterations for SIC
rho_sic = 0.8;   % proportion of identified users for SIC
% mode selection and enable 
fig_en = 0;      % figure enable (locations of UEs and APs)
pilot_sel = 1;   % 0 - identical pilot;  1 - DCS pilot design

% parameter initialization
cell_num = 3*layer*(layer-1) + 1; % number of cells
K = Kc*cell_num;
Ka = K*Sa;
Tmax = K;
Sa_ce = ang_spread/180;   
delay_max = N_CP/Bw;
Ps_max = 10^(Ps_dBm/10)*1e-3;
Pn = 10^(Psd_n/10)*1e-3*Bw;
A_R = dftmtx(N_BS)/sqrt(N_BS);


%% Simulation Setup
Nsim = 1;
T_set = 10:5:80;


%% Allocate Memory
% performance without SIC
Pe = zeros(numel(T_set), 1);
FA = zeros(numel(T_set), 1);
PM = zeros(numel(T_set), 1);
NMSE = zeros(numel(T_set), 1);
NMSE_ang = zeros(numel(T_set), 1);
% performance with SIC
Pe_sic = zeros(numel(T_set), 1);
FA_sic = zeros(numel(T_set), 1);
PM_sic = zeros(numel(T_set), 1);
NMSE_ang_sic = zeros(numel(T_set), 1);


%% Simulation Runs
for sim = 1:Nsim
    %% Mento Carlo Samples
    % generate the locations of UEs and APs
    d = generate_cell_user(r_max, r_min, layer, Kc, fig_en);
    % path loss
    loss = 10.^( - (loss_a + loss_b.*log10(d))/10);
    % transmit power and activity indicator
    Ps = zeros(cell_num, K);
    act_flag = zeros(K,1);
    for nn = 1:cell_num
        temp = Ps_max/r_max^p_ad.*d(nn, (nn-1)*Kc+1:nn*Kc).^p_ad;  
        Ps(:, (nn-1)*Kc+1:nn*Kc) = temp(ones(cell_num,1), :);
        act_flag((nn-1)*Kc + randperm(Kc,Kc*Sa)) = 1;
    end
    id_act = find(act_flag == 1);
    % pilot signal
    if pilot_sel == 0
       S = sqrt(0.5).*(randn(Tmax,K) + 1i*randn(Tmax,K));
       S = S(:,:,ones(1,P));
    else 
       S = sqrt(0.5).*(randn(Tmax,K,P) + 1i*randn(Tmax,K,P));
    end 
    % generate channel model  variance: Lp
    H = zeros(K, N_BS*cell_num, P);
    H_ang = zeros(K, N_BS*cell_num, P);
    varH = zeros(K, cell_num*N_BS, P);
    for nn = 1:cell_num
        for kk = 1:Ka
            Lp = randi([Lp_min, Lp_max], 1, 1);
            [H_tmp, H_ang_tmp] = generate_mimo_channel_off_grid (N_BS, N_UE, Lp, ang_spread,...
                                 delay_max, var_gain, rho_rdd_r, rho_rdd_t, P, N, Bw);
            H(id_act(kk), (nn-1)*N_BS+1:nn*N_BS, :) = sqrt(Ps(nn,id_act(kk)).*loss(nn,id_act(kk)))...
                                                      .* permute(H_tmp, [2 1 3]);
            H_ang(id_act(kk), (nn-1)*N_BS+1:nn*N_BS, :) = sqrt(Ps(nn,id_act(kk)).*loss(nn,id_act(kk)))...
                                                          .* permute(H_ang_tmp, [2 1 3]);
        end
        tmp = (Lp_min+Lp_max)/2*Ps(nn,:).'.*loss(nn,:).';
        varH(:,(nn-1)*N_BS+1:nn*N_BS, :) = tmp(:, ones(1,N_BS), ones(1,P));
    end
%     figure
%     mesh(abs(H))
%     figure
%     mesh(abs(H_ang))

    % received signal
    Y = zeros(Tmax, N_BS*cell_num, P);
    R = zeros(Tmax, N_BS*cell_num, P);
    for pp = 1:P
        noise = sqrt(Pn/2).*(randn(Tmax,N_BS*cell_num) + 1i.*randn(Tmax,N_BS*cell_num));
        Y(:,:,pp) = S(:,:,pp)*H(:,:,pp) + noise;   
        R(:,:,pp) = S(:,:,pp)*H_ang(:,:,pp) + noise;  
    end
    % power amplification
    Y = PA.*Y;
    R = PA.*R;
    H = PA.*H;
    H_ang = PA.*H_ang;
    varH = PA^2.*varH;
    varN = PA^2.*Pn;
    
    %% SIC-based Active User Detection and Channel Estimation
    for T_ii = 1:numel(T_set)
        T = T_set(T_ii);
        
        % SIC initialization
        aus_sic = [];                        % identified active user set for SIC
        aus_sic_pre = [];
        Y_re = zeros(T, cell_num*N_BS, P);   % residual received signal
        Y_err = zeros(T, cell_num*N_BS, P);  % received signal error
        Y_err_nmse_pre = 0;
        H_est = zeros(K, cell_num*N_BS, P);
        H_ang_est = zeros(K, cell_num*N_BS, P);
        
        % SIC iterations
        for iter = 1:niter_sic
            % SIC: Remove the signals from the identified users 
            aus_sic = aus_sic(1: ceil(rho_sic*numel(aus_sic)));
            for pp = 1:P      
                for nn = 1:cell_num
                    H_est(:, (nn-1)*N_BS+1:nn*N_BS, pp) = H_ang_est(:, (nn-1)*N_BS+1:nn*N_BS, pp)*A_R.';
                end
                Y_re(:,:,pp) = Y(1:T,:,pp) - S(1:T,aus_sic,pp)*H_est(aus_sic,:,pp);
                Y_err(:,:,pp) = Y(1:T,:,pp) - S(1:T,:,pp)*H_est(:,:,pp);
            end 
            H_re = H;
            H_re(aus_sic,:,:) = H(aus_sic,:,:) - H_est(aus_sic,:,:);
            Y_err_nmse = norm(Y_err(1:end))^2 / sum(sum(sum(abs(Y(1:T,:,:)).^2)));
            
            % spatial domian based active user detection     
            Sa_aud = (Ka-numel(aus_sic))/K;
            [H_amp, lambda, iter_aud, v] = gmmv_amp(Y_re, S(1:T,:,:), damp_aud, niter_amp,...
                         tol_iter, tol_re, nns_sel_aud, Nco, cell_num, Sa_aud, varH, varN, d, H_re); 
%             % for debug
%             H_reduce = zeros(K, cell_num*N_BS, P);
%             H_reduce(aus_sic,:, :) = H_est(aus_sic,:,:);
%             figure
%             mesh(abs(H_reduce))
%             figure
%             mesh(abs(H_re))
%             figure
%             mesh(abs(H_amp))

            % BI-AD
            act_bi = zeros(K,1);
            act_sic = zeros(K,1);
            bi_avg = zeros(K,1);
            for c = 1:cell_num
                id_ue  = (c-1)*Kc+1:c*Kc;
                id_antenna  = (c-1)*N_BS+1:c*N_BS;
                 bi_avg(id_ue) = sum(sum(lambda(id_ue, id_antenna, :), 2), 3) / (N_BS*P);
            end 
%             [bi_avg_sort, id_sort] = sort(bi_avg, 'descend');
            act_bi(find(bi_avg >= th_bi)) = 1;
            act_sic(find(bi_avg >= th_rel)) = 1;
            
            % add identified users
            act_bi(aus_sic_pre) = 1;
            act_sic(aus_sic_pre) = 1;
            % reliable AUS    
            aus_sic = find(act_sic == 1);
            aus_sic_pre = aus_sic;

            % angular domian based channel estiamtion  
            aus_est = find(act_bi == 1);
            [H_ang_amp, ~, iter_ce] = gmmv_amp(R(1:T,:,:), S(1:T,aus_est,:), damp_ce,...
                 niter_amp, tol_iter, tol_re, nns_sel_ce, Nco, cell_num, Sa_ce, varH(aus_est,:,:),...
                 varN, d, H_ang(aus_est,:,:));    
            H_ang_est = zeros(K, cell_num*N_BS, P);
            H_ang_est(aus_est,:,:) = H_ang_amp;
            
            % detection error
            Pe_tmp = sum(abs(act_bi-act_flag)) / K;
            PM_tmp = sum(~act_bi.*act_flag) / Ka;
            FA_tmp = sum(act_bi.*~act_flag) / (K-Ka);
            % NMSE of the spatial domain and angular domain CE
            H_amp(find(act_bi == 0),:,:) = 0;
            err = H_amp - H_re;
            err2 = 0;
            H_norm2 = 0;
            err_ang = H_ang_est - H_ang;
            err2_ang = 0;
            H_ang_norm2 = 0;
            for nn = 1:cell_num
                tmp = err((nn-1)*Kc+1:nn*Kc, (nn-1)*N_BS+1:nn*N_BS, :);
                err2 = err2 + norm(tmp(:))^2;
                tmp = H((nn-1)*Kc+1:nn*Kc, (nn-1)*N_BS+1:nn*N_BS, :);
                H_norm2 = H_norm2 + norm(tmp(:))^2; 
                tmp = err_ang((nn-1)*Kc+1:nn*Kc, (nn-1)*N_BS+1:nn*N_BS, :);
                err2_ang = err2_ang + norm(tmp(:))^2;
                tmp = H_ang((nn-1)*Kc+1:nn*Kc, (nn-1)*N_BS+1:nn*N_BS, :);
                H_ang_norm2 = H_ang_norm2 + norm(tmp(:))^2;
            end
            NMSE_tmp = err2/H_norm2;
            NMSE_ang_tmp = err2_ang/H_ang_norm2;
            
            % compute received signal error
            if iter == 1
               for pp = 1:P
                   Y_err(:,:,pp) = Y(1:T,:,pp) - S(1:T,:,pp)*H_amp(:,:,pp);
               end
               Y_err_nmse = norm(Y_err(1:end))^2 / sum(sum(sum(abs(Y(1:T,:,:)).^2)));
            end
            
            % print status of SIC
            fprintf('SIC processing...\n')
            fprintf('sim = %d, T = %d, iter_sic = %d\n', sim, T, iter);
            fprintf('iter_aud = %d, iter_ce= %d, FA = %4.5f, PM = %4.5f, Pe = %4.5f, NMSE = %7.9f, NMSE_ang = %7.9f, Y_err = %7.9f\n', ...
                     iter_aud, iter_ce, FA_tmp, PM_tmp, Pe_tmp, NMSE_tmp, NMSE_ang_tmp, Y_err_nmse);
             
            % record data without SIC
            if iter == 1
               Pe(T_ii) = Pe(T_ii) + Pe_tmp;
               PM(T_ii) = PM(T_ii) + PM_tmp;
               FA(T_ii) = FA(T_ii) + FA_tmp;
               NMSE(T_ii) = NMSE(T_ii) + NMSE_tmp;
               NMSE_without_sic = NMSE_tmp;
               NMSE_ang(T_ii) = NMSE_ang(T_ii) + NMSE_ang_tmp;
               NMSE_ang_without_sic = NMSE_ang_tmp;
               Pe_without_sic = Pe_tmp;
            end
       
        end

        % record data with SIC
        Pe_sic(T_ii) = Pe_sic(T_ii) + Pe_tmp;
        PM_sic(T_ii) = PM_sic(T_ii) + PM_tmp;
        FA_sic(T_ii) = FA_sic(T_ii) + FA_tmp;
        NMSE_ang_sic(T_ii) = NMSE_ang_sic(T_ii) + NMSE_ang_tmp;    
        
        % print status of SIC
        fprintf('SIC finished...\n')
        fprintf('sim = %d, T = %d, iter_sic = %d\n', sim, T, iter);
        fprintf('FA = %4.5f, PM = %4.5f, Pe = %4.5f, Pe_sic = %4.5f, NMSE = %7.9f, NMSE_ang = %7.9f, NMSE_ang_sic = %7.9f\n', ...
                 FA_tmp, PM_tmp, Pe_without_sic, Pe_tmp, NMSE_without_sic, NMSE_ang_without_sic, NMSE_ang_tmp);
        fprintf('Pe_avg = %4.5f, Pe_sic_avg = %4.5f, NMSE_avg = %7.9f, NMSE_ang_avg = %7.9f, NMSE_ang_sic_avg = %7.9f\n', ...
                 Pe(T_ii)/Nsim, Pe_sic(T_ii)/Nsim, NMSE(T_ii)/Nsim, NMSE_ang(T_ii)/Nsim, NMSE_ang_sic(T_ii)/Nsim);
    end
toc
end


%% Save Data
Para.Kc = Kc;
Para.Sa = Sa;
Para.N_BS = N_BS;
Para.layer = layer;
Para.r_max = r_max;
Para.p = p_ad;
Para.Ps_dBm = Ps_dBm;
Para.P = P;
Para.Lp_min = Lp_min;
Para.Lp_max = Lp_max;
Para.ang_spread = ang_spread;
Para.Nco = Nco;
Para.iter = niter_amp;
Para.damp_aud = damp_aud;
Para.damp_ce = damp_ce;
Para.nns_sel_aud = nns_sel_aud;
Para.nns_sel_ce = nns_sel_ce;
Para.tol_iter = tol_iter;
Para.tol_re = tol_re;
Para.th_bi = th_bi;
Para.th_rel = th_rel;
Para.niter_sic = niter_sic;
Para.rho_sic = rho_sic;
Para.pilot_sel  = pilot_sel;

Perform.Pe = Pe./Nsim;
Perform.PM = PM./Nsim;
Perform.FA = FA./Nsim;
Perform.NMSE = 10.*log10(NMSE./Nsim);
Perform.NMSE_ang = 10.*log10(NMSE_ang./Nsim);

Perform.Pe_sic = Pe_sic./Nsim;
Perform.PM_sic = PM_sic./Nsim;
Perform.FA_sic = FA_sic./Nsim;
Perform.NMSE_ang_sic = 10.*log10(NMSE_ang_sic./Nsim);
save Data_Cloud_Antenna16_P1.mat Nsim T_set Para Perform


%% Figure
% detection error
figure 
semilogy(T_set,Perform.Pe,'k-*', T_set,Perform.PM,'k--x', T_set,Perform.FA,'k-.+', 'linewidth', 1.2);
hold on
semilogy(T_set,Perform.Pe_sic,'r-*', T_set,Perform.PM_sic,'r--x', T_set,Perform.FA_sic,'r-.+', 'linewidth', 1.2);
set(gca, 'GridLineStyle', '-.');
grid on;
xlabel('Time Overhead $G$'); ylabel('\itP_e')

% NMSE
figure 
plot(T_set,Perform.NMSE,'k--*',  T_set,Perform.NMSE_ang,'r-o', T_set,Perform.NMSE_ang_sic,'g-.+', 'linewidth', 1.2);
hold on
set(gca, 'GridLineStyle', '-.');
grid on;
xlabel('Time Overhead $G$'); ylabel('NMSE')





























