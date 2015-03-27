function [ phRecoveredSyms, phRecoveredSyms_v2, phRecoveredSyms_v3, preCrossCorr ] = ...
    qpskSyncBurst( x, Fs, alpha, debugFilename )

phRecoveredSyms_v2 = [];

[preBits_1sps, preSyms_x1] = genPreamble();
preSyms_x2 = interp(preSyms_x1,2);

[eqBurst] = qpskBurstCFOCorrect(x, Fs, debugFilename);   
preCrossCorr_cmplx = xcorr(preSyms_x2',eqBurst');
preCrossCorr = abs(preCrossCorr_cmplx);

[maxVal,maxIdx] = max(preCrossCorr);
preambleIdxStart = length(eqBurst) - maxIdx + 1;
if(preambleIdxStart<1)
    phRecoveredSyms = [];
    preCrossCorr = [];
    return
end
eqIn = eqBurst(preambleIdxStart:end);
eqIn_div2 = eqIn(1:2:end);
[wOpt]= weiner_filter_equalize( eqIn_div2, debugFilename );
wOpt = wOpt./max(abs(wOpt));
wOpt = wOpt.';      % row vectors are nice

%whFilt_unscaled = filter(wOpt,1,eqIn_div2);
% whFilt_unscaled = conv(eqIn_div2, wOpt);
% fft based convolution
fftSize = length(eqIn_div2)+length(wOpt)-1;
eqIn_div2_ext = [eqIn_div2 zeros(1,fftSize-length(eqIn_div2))];
wOpt_ext = [wOpt zeros(1,fftSize-length(wOpt))];
whFilt_unscaled = ifft(fft(eqIn_div2_ext).*fft(wOpt_ext));

whFilt_mean = mean(abs(whFilt_unscaled));
whFilt = whFilt_unscaled./whFilt_mean;
phRecoveredSyms = qpskFirstOrderPLL(whFilt, alpha);
phRecoveredSyms_v2 = qpskFirstOrderPLL_v2(whFilt, alpha);
phRecoveredSyms_v3 = qpskFirstOrderPLL_v3(whFilt, alpha);

if(debugFilename~=0)
    fname = '/tmp/ml_preCrossCorr.txt';
    dlmwrite(fname, preCrossCorr, 'delimiter', '\n');
    
    fname = '/tmp/ml_preCrossCorrIdx.txt';
    dlmwrite(fname, [preambleIdxStart-1], 'delimiter', '\n');       % minus 1 b/c C is 0 based indexing
    
    fname = '/tmp/ml_eqInDecimated.txt';
    dlmwrite(fname, prepareCmplxVecForWrite(eqIn_div2), 'delimiter', '\n');
    
    fname = '/tmp/ml_wOpt_scaled.txt';
    dlmwrite(fname, prepareCmplxVecForWrite(wOpt), 'delimiter', '\n');
    
    fname = '/tmp/ml_whFilt.txt';
    dlmwrite(fname, prepareCmplxVecForWrite(whFilt_unscaled), 'delimiter', '\n');
    
    fname = '/tmp/ml_whFilt_scaled.txt';
    dlmwrite(fname, prepareCmplxVecForWrite(whFilt), 'delimiter', '\n');
    
    fname = '/tmp/ml_phRecoveredSyms.txt';
    dlmwrite(fname, prepareCmplxVecForWrite(phRecoveredSyms_v2), 'delimiter', '\n');
    
end

% if(debugFilename~=0)
%     save(debugFilename, 'eqBurst', 'preCrossCorr', 'eqIn', 'wOpt', ...
%         'whFilt', 'phRecoveredSyms');
% end

end

