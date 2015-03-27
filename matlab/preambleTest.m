% preamble equalization test
[~, d_n] = genPreamble();       % generate preamble @ 1sps
% simulate a channel
h = [1+j 0 3 j];
r = filter(h,1,d_n);

[wfEqOut, w] = weiner_filter_equalize(r);

% lms
r2 = interp(r,2);
lmsOut = equalize_x2(r2);

subplot(2,2,1)
scatter(real(d_n),imag(d_n))
title('Generated Preamble')
subplot(2,2,2)
scatter(real(r),imag(r))
title('Received Preamble Sequence')
subplot(2,2,3)
scatter(real(wfEqOut),imag(wfEqOut))
title('Wiener Filter Equalized')
subplot(2,2,4)
scatter(real(lmsOut),imag(lmsOut))
title('LMS Equalized')