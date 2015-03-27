close all;
clear all;

%Data Generator 
% Define the M-ary number. 
M = 16; 
% There are k bits in a single M-ary QAM symbol. 
k = log2(M); 
 
% Define the symbol frequency Fd and the sampling frequency Fs. 
Fd = 1; 
Fs = 16; 
 
% Give preamble symbols. 
preamble = ([2 6 2 6 2 6 2 6 2 6]); 
 
% Define the number of user data symbols. 
data_len = 10000; 
% Generate a sequence of uniformly distributed random integers in 
% range [0, M-1]. 
data = randint(data_len,1,M); 
 
% Generate a burst packet. 
burst = ([preamble data']); 
 
 
% 16-QAM Mapper 
% Map integers into complex symbols according to the 16-QAM 
% constellation and Gray-encoded in the meantime. 
burst_map = dmodce(burst,Fd,Fd,'qam',M); 
 
 
%Transmit Filter 
% Design an FIR RRC filter. 
% The roll-off factor of the RRC filter is 0.5. 

rolloff = 0.5; 
% The default time delay of the RRC filter is 3. 
delay = 3; 
% The input digital signal has frequency Fd. the sampling frequency 
% for the filter is Fs. 
[num den] = rcosine(Fd, Fs, 'fir/sqrt', rolloff, delay); 
 
% Filter the signal to obtain the transmit waveforms. 
[I_TxWaveform, t] = rcosflt(burst_map(:,1), Fd, Fs, 'filter', num); 
[Q_TxWaveform, t] = rcosflt(burst_map(:,2), Fd, Fs, 'filter', num); 
 
% Combine the I and Q components into a complex signal. 
s_TxWaveform = I_TxWaveform + j*Q_TxWaveform; 
 
 
%AWGN Channel 
% Generate AWGN. 
% Assume the average SNR per bit EbNodB (in dB) is 10dB. 
EbNodB = 10; 
EbNo = 10^(EbNodB/10); 
 
% Convert EbNo to EsNo. 
EsNo = EbNo*k; 
 
% Define the average power of 16-QAM symbols. 
Es = 10; 
 
% The noise variance. 
noise_var = Es/EsNo*(1/2); 
% The standard deviation of the noise. 
noise_std = sqrt(noise_var); 
 
% Model I and Q components of AWGN with zero mean and No/2 variance. 
noise_I = randn(length(s_TxWaveform),1)*noise_std; 
noise_Q = randn(length(s_TxWaveform),1)*noise_std; 
 
% The baseband complex AWGN. 
noise = noise_I + j*noise_Q; 
 
% Add noise to the transmitted signal. 
r_RxWaveform = s_TxWaveform + noise; 
 
 
%Receive Filter 
% Filter the received signal. 
r_FilterOutput = filter(num, den, r_RxWaveform); 
 
 
%DA Frequency Estimator 
% Define the number of samples used for L_F and N_F. 
L_F = 50*Fs; 
N_F = 51*Fs; 
 
% Calculate time average of the filter output. 
r_aver = sum(r_FilterOutput(L_F+1:L_F+N_F).*conj(r_FilterOutput(1:N_F)))/N_F; 
 
% Frequency offset estimation. 
omega = (atan2(imag(r_aver),real(r_aver)))/L_F; 
 
% Phase derotation. 
n = L_F+N_F+1:length(r_FilterOutput); 
x = r_FilterOutput(n).*(exp(-j*omega *n).'); 
 
 
%DA Phase Estimator 
% Define the number of samples used for N_P. 
N_P = 10*Fs; 
 
% Calculate time average of frequency-adjusted signal x. 
x_aver = sum(x(1:N_P))/N_P; 
 
% Phase offset estimation. 
phi = (atan2(imag(x_aver),real(x_aver))); 
 
% Phase derotation. 
n = N_P+1:length(x); 
y = x(n)*exp(-j*phi); 
 
 
%NDA Phase Tracker 
% Convert the sampling rate to the symbol rate. 
y_symbols = y(Fs:Fs:length(y)); 
 
% Define the values of P, M, and N_W. 
L = 16; 
M = 16; 
N_W = 97; 
 
% NDA phase estimation. 
% Calculate the magnitude (ρ ) and the phase (phi) of y_symbols. 
rho = abs(y_symbols(1:N_W)); 
phi = atan2(imag(y_symbols(1:N_W)),real(y_symbols(1:N_W))); 
 
% Apply a nonlinearity factor L to ρ and a multiplication factor M 
% to Φ. 
y_track = rho .^L.*exp(j*M*phi); 
 
% Phase estimation. 
phi = (1/M)*atan2(sum(imag(y_track)),sum(real(y_track))); 
 
% Phase derotation. 
z = y_track*exp(-j*phi ); 

Demapper and Detector 
data_recov = ddemodce(z,Fd,Fd,'qam',M); 
 
 
BER Counter 
[BER_NUM BER_RATIO] = biterr(data,data_recov,k);
