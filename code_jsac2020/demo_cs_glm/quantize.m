function sig_quantize = quantize(sig, N_bits)
% uniform scalar quantization
[m,n] = size(sig);
sig = sig(:);

delta = 60/(2^N_bits);  % quantization interval

b = (-2^N_bits/2 + 1):2^N_bits/2;
codebook = (-1/2 + b).*delta;

index = 2^N_bits/2 + sign(sig).*min(2^N_bits/2,ceil(abs(sig)/delta)) + ceil((1-sign(sig))/2);
sig_quantize = reshape(codebook(index),m,n);
end

