function [H_f, H_ang_f] = generate_FSF_channel_grid_match(N_BS, N, P, Bw, fs, L_cp, sigma2_alpha, Lp)
% generate the uplink frequency selective fading channel, the BS is
% equipped with ULA while single-antenna is considered at user device.

% written by Malong Ke (kemalong@bit.edu.cn), Beijing Institute of Technology
% version: 2019.03.04

% Inputs��
%   N_BS��number of BS antennas
%   N: number of subcarriers
%   P: pilot subcarriers
%   Bw: bandwidth
%   fs: frequency of sampling
%   fc: frequency of carrier
%   L_cp: length of cyclic prefix
%   sigma2_alpha: average power of path, generally, sigma2_alpha = 1;
%   Lp��number of path

% Outputs��
%   H_f��frequency domain channel matrix for pilots, N_BS��P
%   H_ang_f��virtual-angular domain channel matrix for pilots, N_BS��P

%% Generate the AoAs and steering vectors
n_BS = (0:N_BS-1).';	
index_BS = randi(N_BS,[1,1]);
index_BS = index_BS:index_BS+Lp-1; 
index_BS(index_BS>N_BS) = index_BS(index_BS>N_BS) - N_BS;
theta_BS = index_BS/N_BS;
A_BS = exp(-1i*2*pi*n_BS*theta_BS);

%% Generate delay and gains of path
tau_max = L_cp/fs;                
tau = sort(tau_max.*rand(Lp,1));         
alpha = sqrt(sigma2_alpha/2).*(randn(Lp,1)+1i.*randn(Lp,1));  
alpha = sort(alpha,'descend');
index = randperm(Lp);
tau = tau(index);
alpha = alpha(index);
tau_temp = tau(:,ones(1,P));
alpha_temp = alpha(:,ones(1,P));

%% virtual-angular domain transformation matrix
A_R = dftmtx(N_BS)/sqrt(N_BS);

%% channel
p = 1:P;
f_temp = -Bw/2 + Bw.*(p.*N./P-1)./N;
f_temp = f_temp(ones(Lp,1),:);
H_f = A_BS*(alpha_temp.*exp(-1i*2*pi.*tau_temp.*f_temp));
H_ang_f = A_R'*H_f;

end