clc; clear; close all

N = 1024;
rho = 0.1;
M = 512;
nvar = 1e-3;
Qbits = 3;



x = zeros(N,1);
index = randperm(N,ceil(N*rho));
x(index) = sqrt(0.5).*(randn(length(index),1) + 1i.*randn(length(index),1));
A = sqrt(0.5).*(randn(M,N) + 1i.*randn(M,N));
y = A*x + sqrt(nvar/2).*(randn(M,1) + 1i.*randn(M,1));

yQr = quantize(real(y), Qbits);
yQi = quantize(imag(y), Qbits);
yq = yQr + 1i.*yQi;
% figure
figure
stem(yQr, 'r-o')
hold on
stem(real(y), 'b*')

[xHat, lambda] = gmmv_amp(yq, A, 0.3, 200, 1e-5, 0, nvar);

norm(x-xHat).^2/norm(x).^2
% figure
figure
stem(real(x), 'r-o')
hold on
stem(real(xHat), 'b*')



% Initialization for nonlinear case
lar_num = 1e6;
sma_num = 1e-6;
z_A_ext = zeros(M,1);
v_A_ext = lar_num;
gamma2k = sma_num;
r2k = zeros(N,1);

wvar_hat = 1;

delta = 60/2^Qbits;

iterTurbo = 10;
%% GLM recovery
for iter = 1:iterTurbo
    % component-wise MMSE
    [z_B_post_r, v_B_post_r] = outputUpdate(real(yq), real(z_A_ext), v_A_ext*ones(M,1), sqrt(wvar_hat), Qbits, delta);
    [z_B_post_i, v_B_post_i] = outputUpdate(imag(yq), imag(z_A_ext), v_A_ext*ones(M,1), sqrt(wvar_hat), Qbits, delta);
    z_B_post = z_B_post_r + 1i.* z_B_post_i;
    v_B_post = mean(v_B_post_r) + mean(v_B_post_i);    
%     v_B_post = mean(v_B_post);  
    sigma2_tilde = v_B_post.*v_A_ext./(v_A_ext-v_B_post); %  
    sigma2_tilde = lar_num*(sigma2_tilde<0)+sigma2_tilde.*(sigma2_tilde>0);
    sigma2_tilde = min(sigma2_tilde,lar_num);
    sigma2_tilde = max(sigma2_tilde,sma_num);
    y_tilde = sigma2_tilde.*(z_B_post./v_B_post-z_A_ext./v_A_ext);  %  
    sigma2_tilde = mean(sigma2_tilde);
    
    
    [x_amp, lambda, v] = gmmv_amp(y_tilde, A, 0.3, 20, 1e-5, 0, sigma2_tilde);
    z_A_post = A*x_amp;
    v_A_post = mean(abs(A).^2*v);
    v_A_ext = v_A_post.*sigma2_tilde./(sigma2_tilde-v_A_post);
    v_A_ext = lar_num*(v_A_ext<0)+v_A_ext*(v_A_ext>0);
    v_A_ext = min(v_A_ext,lar_num);
    v_A_ext = max(v_A_ext,sma_num);
    z_A_ext = v_A_ext.*(z_A_post./v_A_post-y_tilde./sigma2_tilde);
end


figure
stem(real(x), 'r-o');
hold on
stem(real(x_amp), 'b-*');

NMSE = norm(x-x_amp)^2/norm(x)^2












