function [ w ] = weiner_filter_equalize( eqInput, debugFilename )

% generate the preamble sequence that will be used for training
% expects input @ 1sps

[~, d_n] = genPreamble();

% % numEqTaps = length(d_n);

x_n = eqInput(1:length(d_n));
x_n = x_n(:);

% generate the input correlation matrix
X = fft(x_n,2^nextpow2(2*size(x_n,1)-1));
X_magSq = abs(X).^2;
rxx_ifft = ifft(X_magSq);
m = length(x_n);
rxx = rxx_ifft./m; % Biased autocorrelation estimate
% creates: http://en.wikipedia.org/wiki/Autocorrelation_matrix
toeplitzMatCol = rxx(1:m);
toeplitzMatRow = conj(rxx(1:m));
R = toeplitz(toeplitzMatCol,toeplitzMatRow);

% make the P vector
xc = xcorr(d_n, x_n);
P_row = xc(1:m);
P = P_row(:);

% solve the optimal weights problem
w = R\P;

if(debugFilename~=0)
    fname = '/tmp/ml_dofFFTOutput.txt';
    dlmwrite(fname, prepareCmplxVecForWrite(X), 'delimiter', '\n');
    
    fname = '/tmp/ml_dofIFFTInput.txt';
    dlmwrite(fname, X_magSq, 'delimiter', '\n');
    
    fname = '/tmp/ml_dofIFFTOutput.txt';
    % scale b/c gnuradio does a ifft scaling factor of the size, and does a
    % conjugate too!! wow.
    dlmwrite(fname, prepareCmplxVecForWrite(conj(rxx_ifft).*128), 'delimiter', '\n');        
    
    fname = '/tmp/ml_dofToeplitzCol.txt';
    dlmwrite(fname, prepareCmplxVecForWrite(toeplitzMatCol), 'delimiter', '\n');
    fname = '/tmp/ml_dofToeplitzRow.txt';
    dlmwrite(fname, prepareCmplxVecForWrite(toeplitzMatRow), 'delimiter', '\n');
    
    fname = '/tmp/ml_dofxc.txt';
    dlmwrite(fname, prepareCmplxVecForWrite(xc), 'delimiter', '\n');
    
    fname = '/tmp/ml_P.txt';
    dlmwrite(fname, prepareCmplxVecForWrite(P_row), 'delimiter', '\n');
    
    fname = '/tmp/ml_wOpt.txt';
    dlmwrite(fname, prepareCmplxVecForWrite(w), 'delimiter', '\n');
end

end
