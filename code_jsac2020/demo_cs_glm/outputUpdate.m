function [z_post, vz_post] = outputUpdate(y, z, mz, sigma, NumBits, delta)
% Performs output node update.
%
% NOTE: This function can potentially run into numerical erros. This is due
% to the sub-function evaluateTotalMoment, which performs integration 
% of a gaussian in some integral given by quantizer boundaries. In case
% when this inteval is far from the mean of the normal and the normal has a
% small variance moments might result in 0, although in reality they
% represent some small values, ratio of which is definetely non-zero.

% length of the signal to estimate
m = size(y, 1);

% Total effective noise (AWGN + estiamtion)
mtv = mz + (sigma^2);

% Initialize outputs

% comupte the lower and up bounds
r_low = y - delta/2;
r_low(r_low < -(2^NumBits-1/2)*delta) = -1e50;

r_up = y + delta/2;
r_up(r_up > (2^NumBits-1/2)*delta) = 1e50;
 
% % complex-valued case
% ita1 = (sign(y).*z - min(abs(r_low),abs(r_up)))./sqrt(2*mtv);
% ita2 = (sign(y).*z - max(abs(r_low),abs(r_up)))./sqrt(2*mtv);
% 
% z_post = z + sign(y).*mz./sqrt(mtv).*((normpdf(ita1) - normpdf(ita2))./(normcdf(ita1) - normcdf(ita2)));
% vz_post = mz/2 - mz.^2./(2*mtv).*((ita1.*normpdf(ita1) - ita2.*normpdf(ita2))./(normcdf(ita1) - normcdf(ita2)) + ((normpdf(ita1) - normpdf(ita2))./(normcdf(ita1) - normcdf(ita2))).^2);

% real-valued case

ita1 = (sign(y).*z - min(abs(r_low),abs(r_up)))./sqrt(mtv/2);
ita2 = (sign(y).*z - max(abs(r_low),abs(r_up)))./sqrt(mtv/2);


A = normpdf(ita1) - normpdf(ita2);
B = normcdf(ita1) - normcdf(ita2);
C = ita1.*normpdf(ita1) - ita2.*normpdf(ita2);


D = A./B;  

E = C./B + (A./B).^2; 


Small_toc = 1e-50;
D(abs(B)<Small_toc) = - ita1(abs(B)<Small_toc);
E(abs(B)<Small_toc) = 1;

z_post = z + sign(y).*mz./sqrt(2*mtv).*D;
vz_post = mz/2 - mz.^2./(2*mtv).*(E);
end