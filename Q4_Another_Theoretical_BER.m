clear all;

SNR = 0:30;
nvar = 0.5.*10.^(-SNR./10);

map_k_1 = genQamMap(1); 
map_k_2 = genQamMap(2); 
map_k_3 = genQamMap(3); 

bits_k_1 = rand(1,2*10^5)<0.5;
bits_k_2 = rand(1,4*10^5)<0.5;
bits_k_3 = rand(1,6*10^5)<0.5;

sym_k_1 = qamMod(bits_k_1,map_k_1);
sym_k_2 = qamMod(bits_k_2,map_k_2);
sym_k_3 = qamMod(bits_k_3,map_k_3);

for t = 1:31;
    Noise_r = normrnd(0,nvar(t)^0.5,1,10^5);
    Noise_i = normrnd(0,nvar(t)^0.5,1,10^5);
    Noise = Noise_r+j*Noise_i;

    r_k_1 = sym_k_1+Noise;
    r_k_2 = sym_k_2+Noise;
    r_k_3 = sym_k_3+Noise;

    llr_k_1 = qamSlice(r_k_1,map_k_1,nvar(t));
    llr_k_2 = qamSlice(r_k_2,map_k_2,nvar(t));
    llr_k_3 = qamSlice(r_k_3,map_k_3,nvar(t));

    bits_r_k_1 = llr_k_1>0;
    bits_r_k_2 = llr_k_2>0;
    bits_r_k_3 = llr_k_3>0;

    bersim_k_1(t) = sum(bits_r_k_1~=bits_k_1)/(2*10^5);
    bersim_k_2(t) = sum(bits_r_k_2~=bits_k_2)/(4*10^5);
    bersim_k_3(t) = sum(bits_r_k_3~=bits_k_3)/(6*10^5);
end

rs = 10.^(SNR./20);
berest_k_1 = qfunc(rs);
berest_k_2 = (1/4).*(3.*qfunc((5^-0.5).*rs)+2.*qfunc(((9/5)^0.5).*rs)-qfunc((5^0.5).*rs));
berest_k_3 = (1/12).*(7.*qfunc((1*(21^-0.5)).*rs)+6.*qfunc((3*(21^-0.5)).*rs)-qfunc((5*(21^-0.5)).*rs)+qfunc((9*(21^-0.5)).*rs)-2.*qfunc((11*(21^-0.5)).*rs)+1.*qfunc((13*(21^-0.5)).*rs));

semilogy(SNR,bersim_k_1,'-',SNR,bersim_k_2,'-',SNR,bersim_k_3,'-');
hold on;
semilogy(SNR,berest_k_1,'o',SNR,berest_k_2,'o',SNR,berest_k_3,'o');
set(gca, 'Fontsize', 16);
xlabel('SNR (dB)');
ylabel('BER');
axis([0 30 1e-4 1]);
hold off;
legend('QPSK(Sim)','16QAM(Sim)','64QAM(Sim)','QPSK(Est)','16QAM(Est)','64QAM(Est)');