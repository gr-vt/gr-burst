% read gr debug data

% fname = '/tmp/gr_inputData.txt';
% inputData = dlmread(fname);
% inputData = inputData(1:2:end) + j*inputData(2:2:end);
% inputData = inputData.';



load good_burst
inputData = x;

Fs = 100e3;
debugFilename = 0;

[pll_v1_out, pll_v2_out, pll_v3_out, preCrossCorr] = ...
        qpskSyncBurst(inputData, Fs, .002, debugFilename);
    
figure;
subplot(2,2,1);
plot(preCrossCorr)
grid on
title('Burst Preamble Cross Correlation')

subplot(2,2,2);
scatter(real(pll_v1_out),imag(pll_v1_out))
grid on
title('Constellation w/ v1 PLL')

subplot(2,2,3);
scatter(real(pll_v2_out),imag(pll_v2_out))
grid on
title('Constellation w/ v2 PLL')

subplot(2,2,4)
scatter(real(pll_v3_out),imag(pll_v3_out))
grid on
title('Constellation w/ v3 PLL')
