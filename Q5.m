clear all;
r = -2:0.1:2;
x = 1:41;

map_1 = genQamMap(1); 
llr_1 = qamSlice(r,map_1,0.5);
llr_1_1 = llr_1(2.*x-1);

map_2 = genQamMap(2); 
llr_2 = qamSlice(r,map_2,0.5*10^-0.5);
llr_2_1 = llr_2(4.*x-3);
llr_2_2 = llr_2(4.*x-2);

map_3 = genQamMap(3); 
llr_3 = qamSlice(r,map_3,0.05);
llr_3_1 = llr_3(6.*x-5);
llr_3_2 = llr_3(6.*x-4);
llr_3_3 = llr_3(6.*x-3);

figure(2);
subplot(331);
plot(r,llr_1_1);
axis([-2 2 -20 20]);
set(gca,'ytick', [0]);
grid on; title('k=1, SNR=0dB');xlabel('r');ylabel('LLR(1)');

subplot(334);
plot(r,llr_2_1);
axis([-2 2 -20 20]);
set(gca,'ytick', [0]);
grid on; title('k=2, SNR=5dB');xlabel('r');ylabel('LLR(1)');

subplot(335);
plot(r,llr_2_2);
axis([-2 2 -20 20]);
set(gca,'ytick', [0]);
grid on; title('k=2, SNR=5dB');xlabel('r');ylabel('LLR(2)');

subplot(337);
plot(r,llr_3_1);
axis([-2 2 -20 20]);
set(gca,'ytick', [0]);
grid on; title('k=3, SNR=10dB');xlabel('r');ylabel('LLR(1)');

subplot(338);
plot(r,llr_3_2);
axis([-2 2 -20 20]);
set(gca,'ytick', [0]);
grid on; title('k=3, SNR=10dB');xlabel('r');ylabel('LLR(2)');

subplot(339);
plot(r,llr_3_3);
axis([-2 2 -20 20]);
set(gca,'ytick', [0]);
grid on; title('k=3, SNR=10dB');xlabel('r');ylabel('LLR(3)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

llr_1_1 = log(10.^(llr_1_1));

llr_2_1 = log(10.^(llr_2_1));
llr_2_2 = log(10.^(llr_2_2));

llr_3_1 = log(10.^(llr_3_1));
llr_3_2 = log(10.^(llr_3_2));
llr_3_3 = log(10.^(llr_3_3));

figure(3);
subplot(331);
plot(r,llr_1_1);
axis([-2 2 -20 20]);
set(gca,'ytick', [0]);
grid on; title('k=1, SNR=0dB');xlabel('r');ylabel('LLR(1)');

subplot(334);
plot(r,llr_2_1);
axis([-2 2 -20 20]);
set(gca,'ytick', [0]);
grid on; title('k=2, SNR=5dB');xlabel('r');ylabel('LLR(1)');

subplot(335);
plot(r,llr_2_2);
axis([-2 2 -20 20]);
set(gca,'ytick', [0]);
grid on; title('k=2, SNR=5dB');xlabel('r');ylabel('LLR(2)');

subplot(337);
plot(r,llr_3_1);
axis([-2 2 -20 20]);
set(gca,'ytick', [0]);
grid on; title('k=3, SNR=10dB');xlabel('r');ylabel('LLR(1)');

subplot(338);
plot(r,llr_3_2);
axis([-2 2 -20 20]);
set(gca,'ytick', [0]);
grid on; title('k=3, SNR=10dB');xlabel('r');ylabel('LLR(2)');

subplot(339);
plot(r,llr_3_3);
axis([-2 2 -20 20]);
set(gca,'ytick', [0]);
grid on; title('k=3, SNR=10dB');xlabel('r');ylabel('LLR(3)');