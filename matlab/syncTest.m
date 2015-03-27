% experiment w/ burst stuff
% load burst_capture      % loads x1,x2,x3,x4,x5,x6
% load burst_capture_with_channel_model
% load burst_capture_with_channel_model_large_burst
% 
% x_in = x6;
% h = [1+j 0 3 j];
% r = filter(h,1,x_in);

% load burst_capture_with_channel_model_lots_of_bursts

Fs = 100e3;             % the sample rate of capture

% r = x5;
% r = awgn(r, 20, 'measured');

fid = fopen('/data/e300collect_100ksps.bin', 'r');
x20 = fread(fid, Inf, 'float32');
if(mod(length(x20),2)==1)
    x20 = x20(1:end-1);
end
x20 = x20(1:2:end) + j*x20(2:2:end);
r = x20(31280:59000).';

x_filt = qpskSyncBurst(r, Fs, .002, 0); 
figure;
subplot(2,1,1); 
scatter(real(r),imag(r)); 
subplot(2,1,2); 
scatter(real(x_filt),imag(x_filt))
grid on