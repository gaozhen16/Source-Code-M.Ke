function [Xhat, lambda, iter, v] = gmmv_amp(Y, Phi, damp, niter, tol, tol_re, nns_sel, Nco,...
                                         cell_num, Sa, varH, varN, d, X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GMMV-AMP algorithm for GMMV CS problem (estimating 3D matrix)
% where incremental EM algorithm is used to learn unknown hyper-parameters.

% Inputs:
%   Y: received signal
%   Phi: measurement matrix
%   damp: damping parameter
%   niter: number of AMP iterations
%   tol&tol_re: termination threshold
%   nns_sel: select nearest neighbor set for NNSPL 0: strctured sparsity   1: clustered sparsity 

% Outputs:
%   Xhat: the estimated matrix
%   lambda: belief indicators
%   iter: number of iterations

% Written by Malong Ke (kemalong@bit.edu.cn), Beijing Institute of Technology
% version: 2019.12.06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
% problem dimension
[M, Q, P] = size(Y);
[~, N, ~] = size(Phi);
N_BS = Q/cell_num;
Kc = N/cell_num;
% sparsity ratio 
lambda = Sa.*ones(N,Q,P);
% mean and variance of traget signal, and noise variance 
xmean = 0;
Xhat = xmean.*ones(N,Q,P);
xvar = varH;
v = xvar;
% nvar = varN;
snr0 = 100; 
nvar = norm(Y(:))^2/(snr0+1)/(M*Q*P);
V = ones(M,Q,P);
Z = Y;
% nearest neighboor sparsity pattern learning
index = ones(N,Q,P);
index_left = [zeros(N,1,P), index(:,1:Q-1,:)];
index_right = [index(:,2:Q,:), zeros(N,1,P)];
index_ahead = cat(3, zeros(N,Q,1), index(:,:,1:P-1));
index_latter = cat(3, index(:,:,2:P), zeros(N,Q,1));
index_3D = index_left + index_right + index_ahead + index_latter; % index_up + index_down + 
clear index
clear index_left
clear index_right
clear index_ahead
clear index_latter
% allocate memory
C = zeros(N,Q,P);
D = zeros(N,Q,P);
L_cal = zeros(N,Q,P);
pai = zeros(N,Q,P);
A = zeros(N,Q,P);
B = zeros(N,Q,P);
Yhat = zeros(M,Q,P);
 

%% AMP iterations
for iter = 1:niter
    Xhat_pre = Xhat;
    V_pre = V;
    
    for p = 1:P
        % factor node update
        V(:,:,p) = damp.*V_pre(:,:,p) + (1-damp).*abs(Phi(:,:,p)).^2*v(:,:,p);
        Z(:,:,p) = damp.*Z(:,:,p) + (1-damp).*(Phi(:,:,p)*Xhat(:,:,p) - ...
                                               V(:,:,p).*(Y(:,:,p)-Z(:,:,p))./(nvar+V_pre(:,:,p)));
    
        % variable node update
        D(:,:,p) = 1 ./ ((abs(Phi(:,:,p)).^2).' * (1./(nvar+V(:,:,p))));
        C(:,:,p) = Xhat(:,:,p) + D(:,:,p).*(Phi(:,:,p)'*((Y(:,:,p)-Z(:,:,p))./(nvar+V(:,:,p))));
   
        % posterior mean and variance
        L_cal(:,:,p) = (1/2).*(log(D(:,:,p)./(D(:,:,p)+xvar(:,:,p))) + abs(C(:,:,p)).^2./D(:,:,p)...
                                               - abs(C(:,:,p)-xmean).^2./(D(:,:,p)+xvar(:,:,p))); 
        pai(:,:,p) = lambda(:,:,p) ./ (lambda(:,:,p) + (1-lambda(:,:,p)).*exp(-L_cal(:,:,p))); 
        A(:,:,p) = (xvar(:,:,p).*C(:,:,p)+xmean.*D(:,:,p)) ./ (D(:,:,p)+xvar(:,:,p));
        B(:,:,p) = (xvar(:,:,p).*D(:,:,p)) ./ (xvar(:,:,p)+D(:,:,p));
        Xhat(:,:,p) = pai(:,:,p).*A(:,:,p);
        v(:,:,p) = pai(:,:,p).*(abs(A(:,:,p)).^2 + B(:,:,p)) - abs(Xhat(:,:,p)).^2;
        
        % reconstruct received signal
        Yhat(:,:,p) = Phi(:,:,p)*Xhat(:,:,p);
    end
    
    % EM-based parameter learning
    nvar_temp = abs(Y-Z).^2./abs(1+V./nvar).^2 + V./(1+V./nvar);
    nvar = sum(nvar_temp(:))/M/Q/P;
   
   % exploiting the structured sparsity for refining the update rule of sparsity ratio
   if nns_sel == 0  % each AP performs AUD and CE locally
      for nn = 1:cell_num
          pai_avg = sum(sum(pai((nn-1)*Kc+1:nn*Kc, (nn-1)*N_BS+1:nn*N_BS, :),2),3) ./ (N_BS*P);
          lambda((nn-1)*Kc+1:nn*Kc, (nn-1)*N_BS+1:nn*N_BS, :) = pai_avg(:, ones(1,N_BS), ones(1,P));
      end
      
   elseif nns_sel == 1  % Cloud computing  
         bi_ap = zeros(N,cell_num); % average at AP
         for c = 1:cell_num                                         
             bi_ap(:,c) = sum(sum(pai(:, (c-1)*N_BS+1:c*N_BS, :),2),3) ./ (N_BS*P);
         end
   
         for nn = 1:N
             d_sum = sum(1./d(:,nn),1);
             bi_tmp = 0;
             for c = 1:cell_num
                 bi_tmp = bi_tmp + 1/(d(c,nn)*d_sum)*bi_ap(nn,c);
             end
             lambda(nn, :, :) = bi_tmp;
         end
   
   elseif nns_sel == 2  %  Edge computing new
         bi_ap = zeros(N,cell_num); % average for each APs
         for c = 1:cell_num                                         
             bi_ap(:,c) = sum(sum(pai(:, (c-1)*N_BS+1:c*N_BS, :),2),3) ./ (N_BS*P);
         end
         
         for nn = 1:N
             [~,sel_cell] = sort(d(:,nn));
             sel_cell = sel_cell(1:Nco);
             d_sum = sum(1./d(sel_cell,nn),1);
             bi_tmp = 0;
             for c = 1:numel(sel_cell)
                 bi_tmp = bi_tmp + 1/(d(sel_cell(c),nn)*d_sum)*bi_ap(nn,sel_cell(c));
             end
              lambda(nn, :, :) = bi_tmp;
         end
         
   elseif nns_sel == 3  %  Edge computing
      for nn = 1:N
          [~,sel_cell] = sort(d(:,nn));
          sel_cell = sel_cell(1:Nco);
          
          pai_temp = 0;
          for c = 1:length(sel_cell)
              pai_avg = sum(pai(nn,(sel_cell(c)-1)*N_BS+1:sel_cell(c)*N_BS))/N_BS;
              pai_temp = pai_temp + pai_avg;
          end
       
          for c = 1:length(sel_cell)
              lambda(nn,(sel_cell(c)-1)*N_BS+1:sel_cell(c)*N_BS) = pai_temp/length(sel_cell);
          end
       end
      
   else  % clustered sparsity 
      pai_left = [zeros(N,1,P), pai(:,1:Q-1,:)];
      pai_right = [pai(:,2:Q,:), zeros(N,1,P)];
      pai_ahead = cat(3, zeros(N,Q,1), pai(:,:,1:P-1));
      pai_latter = cat(3, pai(:,:,2:P), zeros(N,Q,1));
      lambda = (pai_left + pai_right + pai_ahead + pai_latter) ./ index_3D;   
   end

   % stopping criteria
   NMSE_iter =  norm(Xhat(:)-Xhat_pre(:))^2 / norm(Xhat_pre(:))^2;
   NMSE_re = norm(Y(:)-Yhat(:))^2 / norm(Y(:))^2;
   if NMSE_iter < tol || NMSE_re < tol_re
       break;
   end
   
end
end







