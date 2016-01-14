function sym = qamMod(bits,map);
%qamMod generate the symbol.
%sym = qamMod(bits,map) returns 1*n row vector sym. (n is described in line 7 below.)
%
%% map is the points in the constellation map for one real dimension of k bits.(It is an 1*2^k row vector.)
%
%% bits is the corresponding bits sequence of the symbol sequence. [It is an 1*2kn row vector. (Every 2k bits corresponding 1 symbol.)]

[m,n] = size(map);
[o,p] = size(bits);
k = log2(n);

if (rem(p/(2*k),1)~=0)
    error('bits should be an 1*2kn row vector. (Where n = 1,2,3... is a positive interger.)');
end

q = 1:p/(2*k);
g = 2*k-1:-1:0;
h = 2.^g;
b = arrayfun(@(x) (h*(bits(2*k*(x-1)+1:2*k*x))'),q);
b = dec2bin(b,2*k);
yr = bin2dec(b(1:p/(2*k),1:k));
yi = bin2dec(b(1:p/(2*k),k+1:2*k));
sym = map(yr+1)+j.*map(yi+1);