% characterizes the sync performance w/ different impediments

fid = fopen('/data/clean_bursts.bin', 'r');
x = fread(fid, Inf, 'float32');
if(mod(length(x),2)==1)
    x = x(1:end-1);
end
x = x(1:2:end) + j*x(2:2:end);
Fs = 100e3;

x = x(59000:76000);
x = x.';
x = x(1001:14700);

evmVec = zeros(1,20);
snrVec = 1:20;
for snr=snrVec
    x_with_noise = awgn(x,snr,'measured');
    x_filt = qpskSyncBurst(x_with_noise, Fs, .002, 0);
    evmVal = qpskEVM(x_filt);
    evmVec(snr) = evmVal;
end

plot(snrVec, evmVec);
grid on
xlabel('SNR')
ylabel('EVM of Synchronizer')
title('Synchronizer Characterization')
