function [ y ] = qpskBurstCFOCorrect( x, Fs, debugMode )

x4 = x.^4;      % take qpsk ^ 4th to extract baudrate lines
% % % X = fftshift(abs(fft(x)));
X4 = fftshift(abs(fft(x4)).^2);     % i square it here b/c it is just a scaling factor, 
                                    % and in C++ its faster to just take real^2+imag^2 
                                    % and not take the square-root
f = linspace(-Fs/2,Fs/2,length(x4));

[maxVal, maxIdx] = max(X4);
cfoEstimate = f(maxIdx)/4;

%%%%% THIS IS ONE WAY TO DO THE CFO CORRECTIOn
to = (0:length(x)-1)/Fs;
freqCorrVector = exp(-j*2*pi*cfoEstimate*to);
y = x.*freqCorrVector;

%%%%% THIS DOES THE EXACT SAME OPERATION, BUT COMPUTATIONALLY BETTER for
%%%%% ARM PROCESSORS
% rdiff = exp(-j*2*pi*cfoEstimate/Fs);
% y = x.*(cumprod(rdiff*ones(1,length(x))));

if(debugMode~=0)
    fname = '/tmp/ml_cfoFFTInput.txt';
    dlmwrite(fname, prepareCmplxVecForWrite(x4), 'delimiter', '\n');
    
    fname = '/tmp/ml_cfoFFTOut_abs.txt';
    dlmwrite(fname, fftshift(X4), 'delimiter', '\n');       % gr doesnt do fftshift, so we have to do this
    
    fname = '/tmp/ml_cfoEstimate.txt';
    dlmwrite(fname, [cfoEstimate], 'delimiter', '\n');
    
    fname = '/tmp/ml_burstCFOCorrected.txt';
    dlmwrite(fname, prepareCmplxVecForWrite(y), 'delimiter', '\n');
end
    
end

