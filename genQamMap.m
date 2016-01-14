function map = genQamMap(nbits); 
%genQamMap generate Qam Map.
%map = genQamMap(nbits) returns 1*2^nbits (1*2^k) row vector map. [2^nbits (2^k) points in the constellation map for one real dimension.]
%
%% nbits is the number of bits for one real dimension. (It is a scalar.)

k = nbits;
A = 1:2:2^k-1;
scale = (2^(k-2)/sum(A.^2))^0.5;
map2 = scale.*(-2^k+1:2:2^k-1);
m=1;
I=0;
for m = 1:k;
    I = [I; bitset(flipud(I),m)];
end
map(I+ones(2^k,1)) = map2;