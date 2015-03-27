%% Process BITS Samples
% x13 - 18 are random message BITS w/preamble
load bits_1bps_txGen
[preBits_1sps, preSyms_1sps] = genPreamble();
figure;
subplot(3,2,1)
plot(abs(xcorr(preBits_1sps,x13)))
title('Bit Correlation -- Burst 1')
grid on
subplot(3,2,2)
plot(abs(xcorr(preBits_1sps,x14)))
title('Bit Correlation -- Burst 2')
grid on
subplot(3,2,3)
plot(abs(xcorr(preBits_1sps,x15)))
title('Bit Correlation -- Burst 3')
grid on
subplot(3,2,4)
plot(abs(xcorr(preBits_1sps,x16)))
title('Bit Correlation -- Burst 4')
grid on
subplot(3,2,5)
plot(abs(xcorr(preBits_1sps,x17)))
title('Bit Correlation -- Burst 5')
grid on
subplot(3,2,6)
plot(abs(xcorr(preBits_1sps,x18)))
title('Bit Correlation -- Burst 6')
grid on

%% Process TX Generated QPSK Samples
% x7 - x12 are QPSK symbols (1 sample per symbol) of generated messages (w/ preamble)
% for efficiency purposes, they are loaded into a mat file
load qpsk_1sps_txGen
[preBits_1sps, preSyms_1sps] = genPreamble();
figure;
subplot(3,2,1)
plot(abs(xcorr(x7,preSyms_1sps)))
title('Symbol Correlation -- Burst 1')
grid on
subplot(3,2,2)
plot(abs(xcorr(x8,preSyms_1sps)))
title('Symbol Correlation -- Burst 2')
grid on
subplot(3,2,3)
plot(abs(xcorr(x9,preSyms_1sps)))
title('Symbol Correlation -- Burst 3')
grid on
subplot(3,2,4)
plot(abs(xcorr(x10,preSyms_1sps)))
title('Symbol Correlation -- Burst 4')
grid on
subplot(3,2,5)
plot(abs(xcorr(x11,preSyms_1sps)))
title('Symbol Correlation -- Burst 5')
grid on
subplot(3,2,6)
plot(abs(xcorr(x12,preSyms_1sps)))
title('Symbol Correlation -- Burst 6')
grid on

%% process the pulse shaped file to see if we can still detect the correlation peaks
% Notice that we didn't pulse shape the preamble to cross-correlate and we
% are still getting good correlation peaks (peaks are twice the strength of
% the rest of the signal
[preBits_1sps, preSyms_1sps] = genPreamble();

fid = fopen('/home/kiran/test2.bin', 'r');
x19 = fread(fid, Inf, 'float32');
x19 = x19(1:2:end) + j*x19(2:2:end);
nonzeroIdx = find(x19~=0);
x19 = x19(nonzeroIdx);

% rrc filter
x19 = rxPulseShape(x19);

preSyms_x2 = interp(preSyms_1sps,2);
figure;
subplot(2,2,[1 2])
plot(abs(x19))
grid on
title('Burst TX GRC Block -- Output of RRC Filter')
subplot(2,2,3)
corrVec = abs(xcorr(x19,preSyms_x2));
title('Correlation with Preamble')
plot(corrVec);
subplot(2,2,4)
scatter(real(x19),imag(x19))
title('IQ Plot of Transmitted Signal')

%% Process the test OTA received collect
clear;
clc;
close all;
fid = fopen('/data/b210_100ksps.bin', 'r');
% fid = fopen('/data/loopback_collect_100ksps.bin', 'r');
% fid = fopen('/data/e300collect3_100ksps.bin', 'r');
% fid = fopen('/data/e300collect4_100ksps.bin', 'r');
x20 = fread(fid, Inf, 'float32');
if(mod(length(x20),2)==1)
    x20 = x20(1:end-1);
end
x20 = x20(1:2:end) + j*x20(2:2:end);
Fs = 100e3;

[preBits_1sps, preSyms_x1] = genPreamble();
preSyms_x2 = interp(preSyms_x1,2);

[burstApproxIdxs] = amplitudeBurstDetector(x20, Fs, 101, 0);
fprintf('Num Bursts Detected = %d\n', length(burstApproxIdxs));
% Process the OTA bursts
h1 = figure;
for ii=1:length(burstApproxIdxs)
    
    mystruct = burstApproxIdxs{ii};
    burst = x20(mystruct.burstStartIdx:mystruct.burstEndIdx);
    x2 = burst.';        % transpose w/out conjugate  

    %     % optimize the PLL parameter
%     bestEVM = inf;
%     bestAlpha = 1;
%     for alpha=0.0001:.0005:.01
%         whFilt2 = qpskSyncBurst(x2, Fs, alpha);
%         evm = qpskEVM(whFilt2);
%         if(evm<bestEVM)
%             bestEVM = evm;
%             bestAlpha = alpha;
%         end
%     end
%     whFilt2 = qpskSyncBurst(x2, Fs, bestAlpha);
%     fprintf('Best Alpha = %f, Associated EVM = %f dB\n', bestAlpha, bestEVM);

%     if(ii<=5)
%         debugFilename = sprintf('/tmp/ml_test%d.mat', ii);
%     else
%         debugFilename = 0;
%     end
    debugFilename = 0;

    [pll_v1_out, pll_v2_out, pll_v3_out, preCrossCorr] = ...
        qpskSyncBurst(x2, Fs, .002, debugFilename);
    
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
        
    pause;
end

close all;