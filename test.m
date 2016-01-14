clear all;
Nrx = 4;
Ntx = 2;
d = 1/2;
xvar = 1; 
nvar = 1;
theta = [20,40];
SNR = [20,5];
phi = 2*pi*rand(1,2);

H = genGainMatrix(Nrx,Ntx,nvar,d,SNR,theta,phi)
[F_mmse,mse_mmse] = linEst(H,xvar,nvar,'mmse')
[F_zf,mse_zf] = linEst(H,xvar,nvar,'zf')
[F_mf,mse_mf] = linEst(H,xvar,nvar,'mf')

Ntest = 2^16;
xr = normrnd(0,xvar/2,Ntx,Ntest);
xi = normrnd(0,xvar/2,Ntx,Ntest);
nr = normrnd(0,nvar/2,Nrx,Ntest);
ni = normrnd(0,nvar/2,Nrx,Ntest);
x = xr+j*xi;
n = nr+j*ni;

y = H*x+n;
xe_mmse = F_mmse*y;
xe_zf = F_zf*y;
xe_mf = F_mf*y;

B = abs(x).^2;
A_mmse = abs(x-xe_mmse).^2;
A_zf = abs(x-xe_zf).^2;
A_mf = abs(x-xe_mf).^2;

for t = 1:Ntx; 
    NMSE_empirical_mmse(t) = 10*log10(sum(A_mmse(t,:))/sum(B(t,:)));
    NMSE_empirical_zf(t) = 10*log10(sum(A_zf(t,:))/sum(B(t,:)));
    NMSE_empirical_mf(t) = 10*log10(sum(A_mf(t,:))/sum(B(t,:)));
end

NMSE_empirical_mmse
NMSE_empirical_zf
NMSE_empirical_mf

NMSE_expexted_mmse = 10*log10(mse_mmse./xvar)
NMSE_expexted_zf = 10*log10(mse_zf./xvar)
NMSE_expexted_mf = 10*log10(mse_mf./xvar)


for t = 0:90;
    theta = [0,t];
    H = genGainMatrix(Nrx,Ntx,nvar,d,SNR,theta,phi);
    [F_mmse,mse_mmse] = linEst(H,xvar,nvar,'mmse');
    [F_zf,mse_zf] = linEst(H,xvar,nvar,'zf');
    [F_mf,mse_mf] = linEst(H,xvar,nvar,'mf');
    NMSE_expexted_mmse_different_theta(t+1,:) = 10*log10(mse_mmse./xvar);
    NMSE_expexted_zf_different_theta(t+1,:) = 10*log10(mse_zf./xvar);
    NMSE_expexted_mf_different_theta(t+1,:) = 10*log10(mse_mf./xvar);
end

t = 0:90;
figure(1);
subplot(231);
plot(t,NMSE_expexted_mmse_different_theta(:,1));title('NMSE(1) versus theta(2) of mmse');xlabel('theta(2) (degree)');ylabel('NMSE(1) (dB)');
subplot(232);
plot(t,NMSE_expexted_zf_different_theta(:,1));;title('NMSE(1) versus theta(2) of zf');xlabel('theta(2) (degree)');ylabel('NMSE(1) (dB)');
subplot(233);
plot(t,NMSE_expexted_mf_different_theta(:,1));;title('NMSE(1) versus theta(2) of mf');xlabel('theta(2) (degree)');ylabel('NMSE(1) (dB)');
subplot(234);
plot(t,NMSE_expexted_mmse_different_theta(:,2));;title('NMSE(2) versus theta(2) of mmse');xlabel('theta(2) (degree)');ylabel('NMSE(2) (dB)');
subplot(235);
plot(t,NMSE_expexted_zf_different_theta(:,2));;title('NMSE(2) versus theta(2) of zf');xlabel('theta(2) (degree)');ylabel('NMSE(2) (dB)');
subplot(236);
plot(t,NMSE_expexted_mf_different_theta(:,2));;title('NMSE(2) versus theta(2) of mf');xlabel('theta(2) (degree)');ylabel('NMSE(2) (dB)');

theta = [20,40];
for t = -10:20;
    SNR = [10,t];
    H = genGainMatrix(Nrx,Ntx,nvar,d,SNR,theta,phi);
    [F_mmse,mse_mmse] = linEst(H,xvar,nvar,'mmse');
    [F_zf,mse_zf] = linEst(H,xvar,nvar,'zf');
    [F_mf,mse_mf] = linEst(H,xvar,nvar,'mf');
    NMSE_expexted_mmse_different_power(t+11,:) = 10*log10(mse_mmse./xvar);
    NMSE_expexted_zf_different_power(t+11,:) = 10*log10(mse_zf./xvar);
    NMSE_expexted_mf_different_power(t+11,:) = 10*log10(mse_mf./xvar);
end

t = -10:20;
figure(2);
subplot(231);
plot(t,NMSE_expexted_mmse_different_power(:,1));title('NMSE(1) versus SNR(2) of mmse');xlabel('SNR(2) (db)');ylabel('NMSE(1) (dB)');
subplot(232);
plot(t,NMSE_expexted_zf_different_power(:,1));title('NMSE(1) versus SNR(2) of zf');xlabel('SNR(2) (db)');ylabel('NMSE(1) (dB)');
subplot(233);
plot(t,NMSE_expexted_mf_different_power(:,1));title('NMSE(1) versus SNR(2) of mf');xlabel('SNR(2) (db)');ylabel('NMSE(1) (dB)');
subplot(234);
plot(t,NMSE_expexted_mmse_different_power(:,2));title('NMSE(2) versus SNR(2) of mmse');xlabel('SNR(2) (db)');ylabel('NMSE(2) (dB)');
subplot(235);
plot(t,NMSE_expexted_zf_different_power(:,2));title('NMSE(2) versus SNR(2) of zf');xlabel('SNR(2) (db)');ylabel('NMSE(2) (dB)');
subplot(236);
plot(t,NMSE_expexted_mf_different_power(:,2));title('NMSE(2) versus SNR(2) of mf');xlabel('SNR(2) (db)');ylabel('NMSE(2) (dB)');