% load good_burst
% load noise_but_good_burst
% load burst_grnotworking
% load burst_grnotworking2
load burst_grnotworking3

Fs = 100e3;
debugFilename = 0;

fname = '/tmp/burst1.txt';
dlmwrite(fname, prepareCmplxVecForWrite(x), 'delimiter', '\n');

[pll_v1_out, pll_v2_out, pll_v3_out, preCrossCorr] = ...
        qpskSyncBurst(x, Fs, .002, debugFilename);
figure;
subplot(2,2,1);
plot(preCrossCorr)
grid on
title('Burst Preamble Cross Correlation')

subplot(2,2,2);
scatter(real(pll_v1_out),imag(pll_v1_out))
grid on
title('Constellation w/ v1 PLL')