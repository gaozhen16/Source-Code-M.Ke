function [Xhat, lambda, iter] = dmmv_amp(Y, S, cell_num, varH, d, lambda, damp, niter, ...
                                nnspl_mod, tol)
% This simulation code package is mainly used to reproduce the results of the following paper [1]:
% [1] M. Ke et al., "Compressive Massive Access for Internet of Things: Cloud Computing or Fog Computing?,"
% in Proc. IEEE International Conference on Communications (ICC), Virtual Conference, Jun. 2020, pp. 1¨C6.

% If you use this simulation code package in any way, please cite the original paper [1] above. 
% The author in charge of this simulation code pacakge is: Malong Ke (email: kemalong@bit.edu.cn).

% Input:
%   Y: received signal
%   S: pilot matrix
%   cell_num: number of samll cells
%   varH: channel variance
%   d: distance between the users and the devices
%   lambda: active ratio
%   damp: damping
%   niter: number of iterations
%   nnspl_mod:
%   tol: stop criterion


%% Initialization
% problem dimension
[L, P] = size(Y);
[~, K] = size(S);
M = P/cell_num;
Kc = K/cell_num;

% hyper-parameter initialization
lambda = lambda.*ones(K,P);
xmean = 0;
Xhat = xmean.*ones(K,P);
xvar = zeros(K,P);
for nn = 1:cell_num
    varH_temp = varH(nn,:).';
    xvar(:, (nn-1)*M+1:nn*M) = varH_temp(:,ones(1,M));
end
v = xvar;
V = ones(L,P);
Z = Y;
snr0 = 100;
nvar = norm(Y, 'fro')^2/(1+snr0)/L/P; 
 

%% main iterations
for iter = 1:niter
    Xhat_pre = Xhat;
    V_pre = V;
    
    % factor node update
    V = damp.*V_pre + (1-damp).*abs(S).^2*v;
    Z = damp.*Z + (1-damp).*(S*Xhat - V.*(Y-Z)./(nvar+V_pre));
    
    % variable node update
    D = 1 ./ ( (abs(S).^2).' * (1./(nvar+V)) );
    C = Xhat + D.*( S' * ((Y-Z)./(nvar+V)) );
   
    % posterior mean and variance
    L_cal = (1/2).*( log(D./(D+xvar)) + abs(C).^2./D - abs(C-xmean).^2./(D+xvar) ); 
    pai = lambda ./ ( lambda + (1-lambda).*exp(-L_cal) ); 
    A = (xvar.*C + xmean.*D) ./ (D+xvar);
    B = (xvar.*D) ./ (xvar+D);
    Xhat = pai.*A;
    v = pai.*(abs(A).^2 + B) - abs(Xhat).^2;   
    
    % noise learning
    nvar_temp = abs(Y-Z).^2./abs(1+V./nvar).^2 + V./(1+V./nvar);
    nvar = sum(nvar_temp(:))/L/P;
   
   % sparsity learning
   if nnspl_mod == 0  % cloud computing
      pai_temp = sum(pai,2)./P;
      lambda = pai_temp(:,ones(1,P));
   else % fog computing             
      Nco = 3;
      for kk = 1:K
          cell_user = ceil(kk/Kc);
          [~,sel_cell] = sort(abs(d(:,kk) - d(cell_user,kk)));
          sel_cell = sel_cell(1:Nco);
          pai_temp = 0;
          for nn = 1:length(sel_cell)
              pai_avg = sum(pai(kk,(sel_cell(nn)-1)*M+1:sel_cell(nn)*M))/M;
              pai_temp = pai_temp + pai_avg;
          end
          for nn = 1:length(sel_cell)
              lambda(kk,(sel_cell(nn)-1)*M+1:sel_cell(nn)*M) = pai_temp/length(sel_cell);
          end
      end
   end

   % stop criterion
   NMSE_stop = norm(Xhat_pre(:)-Xhat(:))^2/norm(Xhat(:))^2;
   if NMSE_stop < tol
      break;
   end 
  
end

end