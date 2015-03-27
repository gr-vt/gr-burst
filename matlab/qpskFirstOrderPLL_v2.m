function [ y ] = qpskFirstOrderPLL_v2( x, alpha )

phiHat = 1;
% alpha = 0.002;
y = zeros(1,length(x));
for ii=1:length(x)
    y(ii) = x(ii)*phiHat;
    
%     % demodulating circuit
%     if(real(y(ii))>=0 && imag(y(ii))>=0)
%         % 1 + 1j;
%         xHat = exp(1j*pi/2);
%     elseif(real(y(ii))>=0 && imag(y(ii))<0)
%         % 1 - 1j;
%         xHat = exp(1j*3*pi/2);
%     elseif(real(y(ii))<0 && imag(y(ii))<0)
%         % -1 - 1j;
%         xHat = exp(1j*5*pi/2);
%     else
%         % -1 + 1j;
%         xHat = exp(1j*7*pi/2);
%     end

    % the 2 demodulating circuits (above and below) are identical, I was
    % just testing stuff when porting to C

    % demodulating circuit
    if(real(y(ii))>=0 && imag(y(ii))>=0)
        % 1 + 1j;
        xHat = sqrt(2)/2+j*sqrt(2)/2;
    elseif(real(y(ii))>=0 && imag(y(ii))<0)
        % -1 + 1j;
        xHat = sqrt(2)/2 - j*sqrt(2)/2;
    elseif(real(y(ii))<0 && imag(y(ii))<0)
        % -1 - 1j;
        xHat = -sqrt(2)/2 - j*sqrt(2)/2;
    else
        % -1 + 1j;
        xHat = -sqrt(2)/2 + j*sqrt(2)/2;
    end
    
    phiHatT = conj(xHat)*y(ii);
    phiHatT = phiHatT/abs(phiHatT);
    phiHat = conj(phiHatT.^alpha) .* phiHat;
end

end

