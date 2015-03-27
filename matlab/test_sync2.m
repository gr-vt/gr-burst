
close all;
clear all;

sym = floor(rand([1 1024]) * 4);


samp = (1:2048);
symi = floor( (samp+1)/2 );
sps = 2;

theta = pi/4;
%theta = 0;

x = exp(j*2*pi*samp/sps + j*theta) .* exp(j*2*pi*sym(symi)/4);
r = awgn(x, 20);




fd = 0.001
r = r .* exp(j*2*pi*fd*(1:length(r)));


%fo = -0.5:0.001:0.5
%fop = fo
%for i = 1:length(fo)
%    fop(i) = sum(abs( r .* conj(exp(j*(2*pi+fo(i))*samp/sps)) - pi));
%end
%
%plot(fo,fop);


%plot(r,'.');
%
%p = r .* exp(j*2*pi*samp/sps - pi );
%th_ml = -atan( sum( imag( p ) / real( p ) ) )

N = 4;
%th_ml = atan( sum( r.^N ) )
Lf = 1;
omega_ml = (1/Lf) * atan( imag( mean( r.*conj(circshift(r,[0 Lf]))) ) / real( mean( r.* conj(circshift(r,[0 Lf])) ) ) )




%x = exp(j*2*pi*  )






