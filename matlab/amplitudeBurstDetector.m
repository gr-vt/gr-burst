% Amplitude burst detector used to kick off processing down the chain
% for the 802.11 detector

function [burstApproxIdxs, inputSigPower, inputSigPowerSmoothed, chunkMeanVec, threshVec, ratioVec, meanChunkBufferSumVec] = ...
    amplitudeBurstDetector(inputSig, inputSigFs, maSize, plotsFlag)
    % inputSig - the signal on which we try to detect bursts
    % maSize - size of the moving average, this must be odd

% first get the energy of the signal, we care about
% energy b/c the signal could be  +/-
inputSigPower = (abs(inputSig).^2)/2;

% perform moving average of the signal
kernel = ones(1,maSize)/maSize;
inputSigPowerSmoothed = conv(inputSigPower, kernel);
inputSigPowerSmoothed = inputSigPowerSmoothed(1:length(inputSigPower));

if plotsFlag
    figure;
    subplot(3,1,1)
    plot(inputSig);
    axis([0 length(inputSig) min(abs(inputSig)) max(abs(inputSig)) ])
    title('Input Signal')
    grid on
    subplot(3,1,2)
    plot(abs(inputSigPower))
    axis([0 length(inputSig) min(inputSigPower) max(inputSigPower)])
    title('Input Signal Energy')
    grid on
    subplot(3,1,3)
    plot(abs(inputSigPowerSmoothed))
    axis([0 length(inputSig) min(inputSigPowerSmoothed) max(inputSigPowerSmoothed)])
    title('Input Signal Energy (Smoothed)')
    grid on
end

% define chunk size, over which we do mean & thresholding
% to detect whether this is a burst or not
% this should depend on the time delta between which we expect bursts.  the
% smaller this gets, the more processor intensive the burst detector gets.
% chunkSizeInTime = 100e-6;
% chunkSizeInSamples = chunkSizeInTime*inputSigFs;
chunkSizeInSamples = 100;

% we do a floor operation here, and inevitably this ignores the
% last chunk which is less than the other chunks in size
numLoopIter = floor(length(inputSigPowerSmoothed)/chunkSizeInSamples);

meanChunkBufferSize = 3;
meanChunkBuffer = zeros(1,meanChunkBufferSize);
alpha = 0.4;    % the higher this value, the more the weighting is on
                % previous meanChunks
dbThresh_high = 3;        % 3dB = twice the power??
dbThresh_low = 0;

burstApproxIdxs = {};
baIdx = 1;

state = 0;
mystruct = struct;
chunkMeanVec = zeros(1,numLoopIter);
threshVec = zeros(1,numLoopIter);
ratioVec = zeros(1,numLoopIter);
meanChunkBufferSumVec = zeros(1,numLoopIter);
for ii=1:numLoopIter
    idxLo = (ii-1)*chunkSizeInSamples+1;
    idxHi = ii*chunkSizeInSamples;
    smoothedDataChunk = inputSigPowerSmoothed(idxLo:idxHi);
    chunkMean = mean(smoothedDataChunk);
    
    % rotate the chunkMeanBuffer and put the latest one in the end spot
    % i.e. idx=1 is the oldest data point, idx=end is the most recent chunk mean
    meanChunkBuffer(1:meanChunkBufferSize-1) = meanChunkBuffer(2:meanChunkBufferSize);
    meanChunkBuffer(meanChunkBufferSize) = chunkMean;
    
    % calculate the thresholding using exponential weighting
    thresh = calcExpAvg(meanChunkBuffer, alpha);
    % the thresh is returned as +inf if the meanChunkBuffer isn't full, this means that
    % one of the limitations of this method is that we will miss a burst if it occurs 
    % before our buffer fills up.  I think this is OK!
    
    % declare burst based on snr
    ratio_inDb = 10*log10(chunkMean/thresh);
    
    
    %%% for debugging purposes
    chunkMeanVec(ii) = chunkMean;
    threshVec(ii) = thresh;
    ratioVec(ii) = ratio_inDb;
    meanChunkBufferSumVec(ii) = sum(smoothedDataChunk);
    
    if(state==0)
        if ratio_inDb > dbThresh_high            
            mystruct.burstStartIdx = idxLo;
            state = 1;
        end
    elseif(state==1)
        if ratio_inDb < dbThresh_low
            mystruct.burstEndIdx = idxLo+200;
            state = 0;
            
            % store it
            burstApproxIdxs{baIdx} = mystruct;
            baIdx = baIdx + 1;
            mystruct = struct;      % clear out the struct
        end
    end
    
end

end

function [expAvg] = calcExpAvg(vec, alpha)

    if(length(find(vec==0))>0)
        expAvg = 2^16;      % detector wont work until buffer is full
    else
        % expects the latest data point to be in the 
        % last index of vec, this is weighted the leeast
        ii = [1:length(vec)];
        expWeightVector(1,ii) = (1-alpha).^(ii-1)*alpha;
        expAvg = sum(expWeightVector.*vec);
        
    end
end