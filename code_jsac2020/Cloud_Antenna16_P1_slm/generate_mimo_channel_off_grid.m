function [H, H_ang] = generate_mimo_channel_off_grid(N_BS, N_UE, Lp, ang_spread, delay_max,...
                                                      var_gain, rho_rdd_r, rho_rdd_t, P, N, Bs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate the uplink broadband MIMO channel, the BS and the UE are
% equipped with a ULA, respectively;
% power leakage caused by off-grid is considered;
% redundant dictionary is used for finely quantized angle;
% one-ring channel model (clustered sparsity);
% P subchannels out of total N subchannels are considered; 

% Inputs£º
%  N_BS£ºnumber of BS antennas
%  N_UE: number of UE antennas
%  Lp£ºnumber of MPCs
%  ang_spread: angular spread in degree
%  delay_amx: maximum delay of MPCs
%  var_gain: average power of MPCs, generally, var_gian = 1
%  rho_rdd_r: redundant dictionary A_R -> rho_rdd*N_BS¡ÁN_BS
%  rho_rdd_t: redundant dictionary A_T -> rho_rdd*N_UE¡ÁN_UE
%  P: number of considered subchannels
%  N: number of subcarriers
%  Bs: bandwidth

% Outputs£º
%  H£ºspace-frequency domain channel N_BS¡ÁB_UE¡ÁP
%  H_ang£ºangular-frequency domain channel rho_rdd*N_BS¡Árho_rdd*B_UE¡ÁP

% written by M.Ke (kemalong@bit.edu.cn), Beijing Institute of Technology
% version: v1-2019.11.20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% generate AOAs, AODs, and the related steering vectors
% cluster center is uniformly selected in [-pi/2, pi/2], AOAs are generated
% within a given angular spread
ang_spread = ang_spread/180*pi;  % degree to rad
cluster_center = unifrnd(-pi/2, pi/2);
AOAs = unifrnd(cluster_center-ang_spread/2, cluster_center+ang_spread/2, 1, Lp);
AOAs(AOAs > pi/2) = AOAs(AOAs > pi/2) - pi;
AOAs(AOAs < -pi/2) = AOAs(AOAs < -pi/2) + pi;
n_BS = (0:N_BS-1).';
A_BS = exp(-1i*pi*n_BS*sin(AOAs));
% AODs uniformly dstributed in [-pi/2, pi/2]
AODs = unifrnd(-pi/2, pi/2, 1, Lp);  
n_UE = (0:N_UE-1).';
A_UE = exp(-1i*pi*n_UE*sin(AODs));

%% generate delay and complex gain of MPCs
delay = unifrnd(0, delay_max, 1, Lp);
gain = sqrt(var_gain/2).*(randn(1,Lp) + 1i.*randn(1,Lp));  

%% space-angular transformation matrix (redundant dictionary)
A_R = dftmtx(rho_rdd_r*N_BS)/sqrt(rho_rdd_r*N_BS);
A_R = A_R(1:N_BS,:);
A_T = dftmtx(rho_rdd_t*N_UE)/sqrt(rho_rdd_t*N_UE);
A_T = A_T(1:N_UE,:);

%% genenrate channel
H = zeros(N_BS,N_UE,P);
H_ang = zeros(rho_rdd_r*N_BS,rho_rdd_t*N_UE,P);  
% Note that when rho_rdd_r/ rho_rdd_t > 1, the result in not really the 
% virtual angular representation of the channels 


for p = 1:P
    fp = -Bs/2 + Bs/N*(p*N/P-1);
    Delta = diag(gain.*exp(-1i*2*pi.*delay.*fp)); 
    H(:,:,p) = A_BS*Delta*A_UE';
    H_ang(:,:,p) = A_R'*H(:,:,p)*A_T; 
end
end