function [ eqOut, errorVec ] = equalize_x2( eqInput )

% generate the preamble sequence that will be used for training

[~, preSyms] = genPreamble();

numEqTaps = 23;
w = zeros(numEqTaps,1);        % initialize taps to zero

x = zeros(numEqTaps,1);
mu = 10;       % the learning rate

% training loop
errorVec = [];
eqOut = [];
jj = 1;
W_mat = [];
kk = 1;
for ii=[1:length(eqInput)]
    rxSym = eqInput(ii);
    % shift
    x(2:end) = x(1:end-1);
    x(1) = rxSym;
    
    y_k = w'*x;
    % lms (update weights)
    if(mod(ii,2)==0)
        if(jj<=length(preSyms))
            % error is between generated preamble chips and received
            err = preSyms(jj) - y_k;
            w = w + 2*mu*conj(err)*x;       % only update weights w/ training symbols
            % for diagnostic purposes, store the error
            errorVec = [errorVec abs(err)];
            W_mat(kk,:) = w;
            kk = kk + 1;
        end

        eqOut = [eqOut y_k];
        jj = jj + 1;
    end
end

% % find the minimum error vector, and use that weight matrix??
% [minErrVal, minErrIdx] = min(errorVec);
% w = W_mat(minErrIdx,:);

% figure;
% plot(errorVec);
% pause;

end

