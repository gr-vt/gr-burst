function [ y ] = qpskBurstPhaseCorrect( x )

phEst = 1/4 * angle( sum( x ).^4 );
phCorrection = phEst;
fprintf('Rotating by %f degrees\n', phCorrection*180/pi);
y = exp(1j*phCorrection)*x;

end

