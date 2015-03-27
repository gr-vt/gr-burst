% generates test data for comparing sync_v4 versus the matlab sync
clear;
clc;

% % fid = fopen('/data/b210_100ksps.bin', 'r');
% fid = fopen('/data/snr_hunt_noise_pay.dat', 'r');
% x20 = fread(fid, Inf, 'float32');
% if(mod(length(x20),2)==1)
%     x20 = x20(1:end-1);
% end
% x20 = x20(1:2:end) + j*x20(2:2:end);
% Fs = 100e3;
% 
% [preBits_1sps, preSyms_x1] = genPreamble();
% preSyms_x2 = interp(preSyms_x1,2);
% 
% [burstApproxIdxs] = amplitudeBurstDetector(x20, Fs, 101, 0);
% 
% mystruct = burstApproxIdxs{4};
% burst1 = x20(mystruct.burstStartIdx:mystruct.burstEndIdx);
% burst1 = burst1.';
% if( mod(length(burst1),2)~=0)
%     burst1 = burst1(1:end-1);
% end
% 

load e300_airburst
Fs = 100e3;

% write out burst 1 so gr can read it and do the same processing
fname = '/tmp/burst1.txt';
dlmwrite(fname, prepareCmplxVecForWrite(burst1), 'delimiter', '\n');

debugFilename = 1;

[pll_v1_out, pll_v2_out, pll_v3_out, preCrossCorr] = ...
        qpskSyncBurst(burst1, Fs, .002, debugFilename);

% fprintf('done\n');
    
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
