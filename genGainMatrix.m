function H = genGainMatrix(Nrx,Ntx,nvar,d,SNR,theta,phi);
%genGainMatrix Generate the gain matrix. 
%H = genGainMatrix(Nrx,Ntx,nvar,d,SNR,theta,phi) returns the Nrx*Ntx GainMatrix H with parameters of Nrx,Ntx,nvar,d,SNR,theta,phi.
%
%% Nrx, Ntx is the number of receivers and transmitters. (They are scalars with the units of count.)
%
%% nvar is the vatiance of the additive noise. (It is a scalar with the same unit as the receive power in linear scale. Generally we can use the unit of power in watt.) 
%  Here I assume that the additive noise is complex Gaussian vector with zero mean vector and variance matrix nvar*I.
%
%% d is the the antenna spacing normalized to wavelengths. (It is a scalar with the unit of 1/wavelength.)
%
%% SNR is the average receive power normalized to nvar in db scale. (It is a 1*Ntx vector with the unit of db.)
%
%% theta is the angle of arrival of the signal from the transmitter. (It is a 1*Ntx vector with the unit of degree.)
%
%% phi is the phase of the wave. (It is a 1*Ntx vector with the unit of rad.)

m = Nrx;
n = Ntx;
H = zeros(m,n);
P = nvar*10.^(SNR/10);

for m = 1:Nrx;
    for n = 1:Ntx;
        H(m,n) = (P(n)^0.5)*exp(2*pi*j*m*d*cosd(theta(n))+j*phi(n));
    end
end