function [F,mse] = linEst(H,xvar,nvar,estType);
%linEst Linear estimator 
%[F,mse] = linEst(H,xvar,nvar,estType) return the linear estimator Ntx*Nrx matrix F and mean square error Ntx*1 vector mse in linear scale
%with the parameters of H, xvar, nvar, estType.
%
%% H is the Gain Matrix. (It is a Nrx*Ntx matrix.)
%
%% xvar and nvar are the variances of  transmitted symbol  and additive noise. 
%  (They are scarlars  with the same unit as the receive power in linear scale. Generally we can use the unit of power in watt.)
%  Here I assume that the transmitted symbol and additive noise are complex Gaussian vector 
%  with zero mean vectors and variance matrix xvar*I and nvar*I respectively.
%
%% The default estType is mmse
%
%% [F,mse] = linEst(H,xvar,nvar) = linEst(H,xvar,nvar,'mmse') return the mmse linear estimator Ntx*Nrx matrix F and mmse mean square error Ntx*1 vector mse in linear scale.
%  [F,mse] = linEst(H,xvar,nvar,'zf') return the zf linear estimator Ntx*Nrx matrix F and zf mean square error Ntx*1 vector mse in linear scale.
%  [F,mse] = linEst(H,xvar,nvar,'mf') return the mf linear estimator Ntx*Nrx matrix F and mf mean square error Ntx*1 vector mse in linear scale.

[m,n] = size(H);
r = xvar/nvar;
e = 10^-6*max(diag(H'*H));
H2 = abs(H').^2;
HC =H';

for x = 1:n;
    for y = 1:n;
        B(x,y) = HC(x,:)*H(:,y);
    end
end

HH = abs(B).^2; 

for t = 1:n;
    a(t) = (xvar*real(HC(t,:)*H(:,t)))/(xvar*sum(HH(t,:))+nvar*sum(H2(t,:)));
end

% Set default estimation type if not specified
if (nargin < 4)
    estType = 'mmse';
end

switch estType
    case 'mmse'
        F = r.*inv(eye(n,n)+r.*H'*H)*H';
    case 'zf'
        F = inv(H'*H+e.*eye(n,n))*H';
    case 'mf'
        F = diag(a)*H';
end

%Compute the mse vector
F2 = F.*conj(F); 

for p = 1:n;
    for k = 1:n;
        FH2(p,k) = abs(F(p,:)*H(:,k))^2;
    end
end

for p = 1:n;
    mse(p) = xvar*(1-2*real(F(p,:)*H(:,p))+sum(FH2(p,:)))+nvar*sum(F2(p,:));
end