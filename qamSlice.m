function llr = qamSlice(r,map,nvar);
%qamSlice generate the bitwise log likelihood ratio (LLR).
%llr = qamSlice(sym,map,nvar) returns 1*2kn row vector llr. (llr(i) is the LLR for bits(i).)
%
%% r is the received symbol, it is the corresponding symbol from transmitted bits add the noise. (It is a 1*n row vector.)
%
%% map is the points in the constellation map for one real dimension of k bits.(It is an 1*2^k row vector.)
%
%% nvar is the variance of the noise (N0/2). (It is a scalar.)

[m,n] = size(map);
r = r.';
[m_s,n_s] = size(r);
sym_r = real(r);
sym_i = imag(r);
sym_r = sym_r*ones(1,n/2);
sym_i = sym_i*ones(1,n/2);
k = log2(n);
y1 = 1:n;
y2 = dec2bin(y1-1)-'0';
y3 = ~y2;
for x = 0:k-1;
    y4 = upsample(1,k,x);
    y5 = diag(y2*y4*y1);
    y6 = diag(y3*y4*y1);
    y5 = y5(y5~=0);
    y6 = y6(y6~=0);   
    map_1 = ones(m_s,1)*map(y5);
    map_0 = ones(m_s,1)*map(y6);
    Z_r_1 = log10(sum(exp((sym_r-map_1).^2./(2*nvar).*(-1)),2));
    Z_r_0 = log10(sum(exp((sym_r-map_0).^2./(2*nvar).*(-1)),2));
    Z_i_1 = log10(sum(exp((sym_i-map_1).^2./(2*nvar).*(-1)),2));
    Z_i_0 = log10(sum(exp((sym_i-map_0).^2./(2*nvar).*(-1)),2));
    llr_r(:,x+1) = Z_r_1-Z_r_0; 
    llr_i(:,x+1) = Z_i_1-Z_i_0;    
end
llr = [llr_r llr_i];
llr = llr';
llr = llr(:);
llr = llr.';