clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
QPSK_Symbol=exp(1i*[pi/4 3*pi/4 5*pi/4 7*pi/4]);                              % Creat QPSK_Symbol 
n=2^16;                                                                       % number of the sequence
x=QPSK_Symbol(randi(4,n,1));                                                  % Creat Sequence x
xUP=upsample(x,2);                                                            % upsample X to Xup
FIR=fir1(60,0.5);                                                             % Creat FIR
xPS=filter(FIR,1,xUP);                                                        % Creat Xps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nfft = 1024;                                                                  % PSD of xPS
fsym = 10e+6;
fsxps = 2*fsym;
[Sps,wps] = pwelch(xPS, [], [], nfft,fsxps);
Sps_db = 10*log10(Sps);
rwps(1:512) = wps(513:1024)-fsxps;
rwps(513:1024) = wps(1:512);
rSps_db(1:512) = Sps_db(513:1024);
rSps_db(513:1024) = Sps_db(1:512);
figure(1);
subplot(221);
plot(rwps,rSps_db);title('PSD of xPS');xlabel('Hz');ylabel('dB');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xDAC=upsample(xPS,4)+upsample(xPS,4,1)+upsample(xPS,4,2)+upsample(xPS,4,3);   % Creat xDAC
ovPS = 2;                                                                     % pulse shape oversampling ratio
ovRF = 4;                                                                     % RF oversampling ratio
nfiltRF = 3;                                                                  % filter order
ripRF = 0.5;                                                                  % passband ripple
wpRF = 1/ovPS/ovRF;                                                           % digital cutoff
[bRF,aRF] = cheby1(nfiltRF,ripRF,wpRF);
xRF = filter(bRF, aRF, xDAC);                                                 % filter the signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                  
fsxdac=4*fsxps;                                                               % PSD of xDAC  
[Sdac,wdac] = pwelch(xDAC, [], [], nfft,fsxdac);
Sdac_db=10*log10(Sdac);
rwdac(1:512) = wdac(513:1024)-fsxdac;
rwdac(513:1024) = wdac(1:512);
rSdac_db(1:512) = Sdac_db(513:1024);
rSdac_db(513:1024) = Sdac_db(1:512);
figure(1);
subplot(222);
plot(rwdac,rSdac_db,'g');title('PSD of xDAC');xlabel('Hz');ylabel('dB');                                                                                               
[Srf,wrf] = pwelch(xRF, [], [], nfft,fsxdac);                                 % PSD of xRF
Srf_db=10*log10(Srf);
rwrf(1:512) = wrf(513:1024)-fsxdac;
rwrf(513:1024) = wrf(1:512);
rSrf_db(1:512) = Srf_db(513:1024);
rSrf_db(513:1024) = Srf_db(1:512);
figure(1);
subplot(224);
plot(rwrf,rSrf_db,'r');title('PSD of xRF');xlabel('Hz');ylabel('dB');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
figure(2);                                                                    % PSD of xDAC and xRF in single plot
plot(rwdac,rSdac_db,'g',rwrf,rSrf_db,'r');
title('PSD of xDAC and xRF');xlabel('Hz');ylabel('dB');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
st = 8;                                                                       % Select a value from 1 to ovPS*ovRF=8
yRX = xRF(st:ovPS*ovRF:end);                                                  % subsample the RF output
figure(1);
subplot(223);
plot(yRX, '.');                                                               % Plot the received samples
title('Received Symbols with best choice of st=8');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3);                                                                    % Plot the received symbols with different st
subplot(241);
st = 1;                                                                       % Select a value from 1 to ovPS*ovRF=8
yRX = xRF(st:ovPS*ovRF:end);                                                  % subsample the RF output
plot(yRX, '.');                                                               % Plot the samples
title('Received Symbols with st=1');
subplot(242);
st = 2;                                                                       % Select a value from 1 to ovPS*ovRF=8
yRX = xRF(st:ovPS*ovRF:end);                                                  % subsample the RF output
plot(yRX, '.');                                                               % Plot the samples
title('Received Symbols with st=2');
subplot(243);
st = 3;                                                                       % Select a value from 1 to ovPS*ovRF=8
yRX = xRF(st:ovPS*ovRF:end);                                                  % subsample the RF output
plot(yRX, '.');                                                               % Plot the samples
title('Received Symbols with st=3');
subplot(244);
st = 4;                                                                       % Select a value from 1 to ovPS*ovRF=8
yRX = xRF(st:ovPS*ovRF:end);                                                  % subsample the RF output
plot(yRX, '.');                                                               % Plot the samples
title('Received Symbols with st=4');
subplot(245);
st = 5;                                                                       % Select a value from 1 to ovPS*ovRF=8
yRX = xRF(st:ovPS*ovRF:end);                                                  % subsample the RF output
plot(yRX, '.');                                                               % Plot the samples
title('Received Symbols with st=5');
subplot(246);
st = 6;                                                                       % Select a value from 1 to ovPS*ovRF=8
yRX = xRF(st:ovPS*ovRF:end);                                                  % subsample the RF output
plot(yRX, '.');                                                               % Plot the samples
title('Received Symbols with st=6');
subplot(247);
st = 7;                                                                       % Select a value from 1 to ovPS*ovRF=8
yRX = xRF(st:ovPS*ovRF:end);                                                  % subsample the RF output
plot(yRX, '.');                                                               % Plot the samples
title('Received Symbols with st=7');
subplot(248);
st = 8;                                                                       % Select a value from 1 to ovPS*ovRF=8
yRX = xRF(st:ovPS*ovRF:end);                                                  % subsample the RF output
plot(yRX, '.');                                                               % Plot the samples
title('Received Symbols with st=8');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bit_rate = log2(4)*fsym                                                       % Caculate the bit rate
M = mean(rSps_db(257:769));                                                   % First method to find fc and percentage                                                 
M_fc = M-40;
fc_positive = interp1(rSps_db(719:819),rwps(719:819),M_fc,'spline');
fc_negative = interp1(rSps_db(207:307),rwps(207:307),M_fc,'spline');
fc_method1 = (fc_positive-fc_negative)/2
occupied_bandwidth_2fc_method1 = 2*fc_method1
percentage_method1 = (2*fc_method1-fsym)/fsym
fc_method2 = 2*fsym*sum(rSps_db>M_fc)/1024/2                                  % Second method to find fc and percentage 
occupied_bandwidth_2fc_method2 = 2*fc_method2
percentage_method2 = (2*fc_method2-fsym)/fsym
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                            diferent filter length                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FIR_fl40=fir1(40,0.5);                                                        % Different filter length = 40 
xPS_fl40=filter(FIR_fl40,1,xUP);        
[Sps_fl40,wps_fl40] = pwelch(xPS_fl40, [], [], nfft,fsxps);
Sps_fl40_db = 10*log10(Sps_fl40);
rwps_fl40(1:512) = wps_fl40(513:1024)-fsxps;
rwps_fl40(513:1024) = wps_fl40(1:512);
rSps_fl40_db(1:512) = Sps_fl40_db(513:1024);
rSps_fl40_db(513:1024) = Sps_fl40_db(1:512);
FIR_fl80=fir1(80,0.5);                                                        % Different filter length = 80 
xPS_fl80=filter(FIR_fl80,1,xUP);        
[Sps_fl80,wps_fl80] = pwelch(xPS_fl80, [], [], nfft,fsxps);
Sps_fl80_db = 10*log10(Sps_fl80);
rwps_fl80(1:512) = wps_fl80(513:1024)-fsxps;
rwps_fl80(513:1024) = wps_fl80(1:512);
rSps_fl80_db(1:512) = Sps_fl80_db(513:1024);
rSps_fl80_db(513:1024) = Sps_fl80_db(1:512);
figure(4);                                                                    % Plot the PSD of xPS with different filter length 
subplot(222);
plot(rwps,rSps_db);
title('PSD of xPS with filter length=60');xlabel('Hz');ylabel('dB');
subplot(221);
plot(rwps_fl40,rSps_fl40_db,'g');
title('PSD of xPS with filter length=40');xlabel('Hz');ylabel('dB');
subplot(223);
plot(rwps_fl80,rSps_fl80_db,'r');
title('PSD of xPS with filter length=80');xlabel('Hz');ylabel('dB');
subplot(224);
plot(rwps_fl40,rSps_fl40_db,'g',rwps,rSps_db,rwps_fl80,rSps_fl80_db,'r');
title('PSD of xPS with filter length=40,60,80');xlabel('Hz');ylabel('dB');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_fl40 = mean(rSps_fl40_db(257:769));                                         % Find fc and percentage of different filter length = 40                                                
M_fc_fl40 = M_fl40-40;
fc_fl40 = 2*fsym*sum(rSps_fl40_db>M_fc_fl40)/1024/2;                                          
occupied_bandwidth_2fc_fl40 = 2*fc_fl40;
percentage_fl40 = (2*fc_fl40-fsym)/fsym
M_fl80 = mean(rSps_fl80_db(257:769));                                         % Find fc and percentage of different filter length = 80                                                
M_fc_fl80 = M_fl80-40;
fc_fl80 = 2*fsym*sum(rSps_fl80_db>M_fc_fl80)/1024/2;                                        
occupied_bandwidth_2fc_fl80 = 2*fc_fl80;
percentage_fl80 = (2*fc_fl80-fsym)/fsym
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 