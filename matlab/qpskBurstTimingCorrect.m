function [ y ] = qpskBurstTimingCorrect( x, inputSps, interpFactor )

xi = interp(x, interpFactor);
samps_per_sym = inputSps*interpFactor;     % 2 is the original samps / sym
% find the best timing lag
% do timing recovery 
bestTau_max = 1;
bestTau_min = 1;
maxSigPow = -inf;
minSigPow = inf;

numTimingEsts = 600;
if(length(xi)<numTimingEsts)
    numTimingEsts = length(xi)-100;
end
bestTauMinVec = zeros(1,numTimingEsts);
bestTauMaxVec = zeros(1,numTimingEsts);
for startIdx=1:numTimingEsts
    for tau=[1:samps_per_sym]
        sig = xi(tau:samps_per_sym:end);
        sigPow = sum(sig.^2);
        if(sigPow > maxSigPow)
            bestTau_max = tau;
            maxSigPow = sigPow;
        end
        if(sigPow < minSigPow)
            bestTau_min = tau;
            minSigPow = sigPow;
        end
    end
    bestTauMinVec(startIdx) = bestTau_min;
    bestTauMaxVec(startIdx) = bestTau_max;
end
%     bestTau = bestTau_max;
bestTau = mod(bestTau_min+samps_per_sym/2, samps_per_sym);
if(bestTau==0)      % edge case w/ matlab indexing
    bestTau = 1;
end

y = xi(bestTau:samps_per_sym:end);


end

