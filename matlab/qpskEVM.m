function [ evm ] = qpskEVM( x )

evmVec = zeros(1,length(x));
for ii=1:length(x)
    sym = x(ii);
    symPow = abs(sym);
    idealSymScaleVal = sqrt(symPow/2);
    if(real(sym)>=0 && imag(sym)>=0)
        idealSym = idealSymScaleVal*(1 + 1j);
        evmVal = abs(sym-idealSym)/abs(idealSym);
        evmVec(ii) = evmVal;
    elseif(real(sym)>=0 && imag(sym)<0)
        idealSym = idealSymScaleVal*(1 - 1j);
        evmVal = abs(sym-idealSym)/abs(idealSym);
        evmVec(ii) = evmVal;
    elseif(real(sym)<0 && imag(sym)<0)
        idealSym = idealSymScaleVal*(-1 -1j);
        evmVal = abs(sym-idealSym)/abs(idealSym);
        evmVec(ii) = evmVal;
    else
        idealSym = idealSymScaleVal*(-1 + 1j);
        evmVal = abs(sym-idealSym)/abs(idealSym);
        evmVec(ii) = evmVal;
    end
end

evm = 10*log10(nanmean(evmVec));

end

