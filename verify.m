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
[F_mmse,mse_mmse] = linEst(H,xvar,nvar,'mmse');
[F_zf,mse_zf] = linEst(H,xvar,nvar,'zf');
[F_mf,mse_mf] = linEst(H,xvar,nvar,'mf');
mse_mmse
mse_zf
mse_mf

A = eye(Ntx)-F_mmse*H;
P = xvar*A*A'+nvar*F_mmse*F_mmse';
mse_mmse_verify = diag(P)
A = eye(Ntx)-F_zf*H;
P = xvar*A*A'+nvar*F_zf*F_zf';
mse_zf_verify = diag(P)
A = eye(Ntx)-F_mf*H;
P = xvar*A*A'+nvar*F_mf*F_mf';
mse_mf_verify = diag(P)