function [ y ] = qpskFirstOrderPLL( x, alpha )

phiHat = 0;
% alpha = 0.002;
y = zeros(1,length(x));
for ii=1:length(x)
    y(ii) = x(ii)*exp(-j*phiHat);
    
    % demodulating circuit
    if(real(y(ii))>=0 && imag(y(ii))>=0)
        % 1 + 1j;
        xHat = exp(1j*pi/2);
    elseif(real(y(ii))>=0 && imag(y(ii))<0)
        % 1 - 1j;
        xHat = exp(1j*3*pi/2);
    elseif(real(y(ii))<0 && imag(y(ii))<0)
        % -1 - 1j;
        xHat = exp(1j*5*pi/2);
    else
        % -1 + 1j;
        xHat = exp(1j*7*pi/2);
    end
    
    phiHatT = angle(conj(xHat)*y(ii));
    phiHat = phiHatT*alpha + phiHat;
    
end

end

