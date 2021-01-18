function [Xhat, lamda, iter, NMSE_watch] = mmv_amp(Y, A, niter, tol)
% The algorithm used to reproduce the results in the paper:
% M. Ke et al., "Compressive massive random access for massive machine-type communications (mMTC)," 
% in Proc IEEE Global Conf. Signal Inform. Process. (GlobalSIP), Anaheim, USA, Nov. 2018, pp. 156-161

snr0 = 100;
damp_AMP = 0.3;
[G, Q, P] = size(Y);
[~, K, ~] = size(A);
del = G/K;

%% parameters initialization
normal_cdf = @(x) 1/2*(1 + erf(x/sqrt(2)));
normal_pdf = @(x) 1/sqrt(2*pi)*exp(-x.^2/2);
alpha_grid = linspace(0,10,1024);
rho_SE = (1 - (2/del)*((1+alpha_grid.^2).*normal_cdf(-alpha_grid)-alpha_grid.*normal_pdf(alpha_grid)))...
         ./(1 + alpha_grid.^2 - 2*((1+alpha_grid.^2).*normal_cdf(-alpha_grid)-alpha_grid.*normal_pdf(alpha_grid)));
rho_SE = max(rho_SE);
lamda = rho_SE*del*ones(K,Q,P); % the initialized soft support for BG prob
xmean = 0;
xvar = 0;
nvar = zeros(1,P);
for p = 1:P
    for q = 1:Q
        nvar(p) = nvar(p) + norm(Y(:,q,p))^2 ./ ((1+snr0)*G);
        xvar = xvar + (norm(Y(:,q,p))^2 - G.*nvar(p)) ./ ((rho_SE*del)*norm(A(:,:,p),'fro')^2);
    end
end
xvar = xvar/Q/P;
nvar = nvar./Q;

Xhat(1:K,1:Q,1:P) = xmean; 
v(1:K,1:Q,1:P) = xvar;
V = ones(G,Q,P);
Z = Y;
S = zeros(K,Q,P);
R = zeros(K,Q,P);
m = zeros(K,Q,P);
L = zeros(K,Q,P);
Vrs = zeros(K,Q,P);
NMSE_watch  = zeros(niter,1);
%% AMP iteration
for iter = 1:niter
    X_pre = Xhat;
    for p = 1:P
        V_new = abs(A(:,:,p)).^2 * v(:,:,p);
        Z_new = A(:,:,p)*Xhat(:,:,p) - (Y(:,:,p) - Z(:,:,p)) ./ (nvar(p) + V(:,:,p)) .* V_new;
        V(:,:,p) = damp_AMP*V(:,:,p) + (1-damp_AMP)*V_new;
        Z(:,:,p) = damp_AMP*Z(:,:,p) + (1-damp_AMP)*Z_new;
        S(:,:,p) = 1 ./ ( (abs(A(:,:,p)).^2).' * (1./(nvar(p) + V(:,:,p))) );
        R(:,:,p) = ( (conj(A(:,:,p))).' * ( (Y(:,:,p) - Z(:,:,p)) ./ (nvar(p) + V(:,:,p)) ) ) .* S(:,:,p) + Xhat(:,:,p);
        L(:,:,p) = (1/2).*(  log(S(:,:,p)./(S(:,:,p)+xvar)) + (abs(R(:,:,p))).^2./(S(:,:,p)) - (abs(R(:,:,p)-xmean)).^2./(S(:,:,p)+xvar)   );
        lamda(:,:,p) = lamda(:,:,p) ./ (lamda(:,:,p)+(1-lamda(:,:,p)).*exp(-L(:,:,p)));

        m(:,:,p) = (xvar.*R(:,:,p) + xmean.*S(:,:,p)) ./ (S(:,:,p) + xvar);
        Vrs(:,:,p) = (xvar.*S(:,:,p)) ./ (xvar + S(:,:,p));
        Xhat(:,:,p) = lamda(:,:,p).*m(:,:,p);
        v(:,:,p) = lamda(:,:,p).*(abs(m(:,:,p)).^2+Vrs(:,:,p)) - (abs(Xhat(:,:,p))).^2;
    
        nvar(p) = sum(sum( abs(Y(:,:,p)-Z(:,:,p)).^2./abs(1+V(:,:,p)./nvar(p)).^2 + V(:,:,p)./ (1 + V(:,:,p)./nvar(p)) ))/(G*Q);
    end
    xmean = sum(sum(sum( lamda.*m ))) / sum(sum(sum( lamda)));
    xvar = (1./sum(sum(sum( lamda )))).* sum(sum(sum( lamda.*(abs(xmean-m).^2+Vrs) )));                 

    %% parameters learning
    if P == 1
       pai = lamda;
       pai_update = sum(pai, 2)./Q; 
       for k = 1:K
            lamda(k,:,:) = pai_update(k);
       end
    else
        pai = permute(lamda, [2 3 1]);
        pai_update = squeeze(sum(sum(pai)))./(Q*P); 
        for k = 1:K
            lamda(k,:,:) = pai_update(k);
        end
    end

    NMSE_watch(iter)  = norm(X_pre(1:end) - Xhat(1:end))/norm(Xhat(1:end));
    if(NMSE_watch(iter) < tol)
       break;
    end   
end
end