function [ preSeqBits, preSyms ] = genPreamble()

preSeqBits = [0 1 0 0 1 0 1 1 1 1 0 1 0 1 1 0 0 1 1 1 1 1 1 0 0 0 0 0 1 0 ...
    0 1 1 0 0 1 1 0 0 1 1 0 1 1 0 0 1 0 0 0 0 0 1 0 0 1 0 1 1 1 1 1 0 1 ...
    0 0 1 1 1 1 0 0 1 1 0 0 0 0 0 0 1 1 1 1 1 1 0 1 0 1 0 1 1 0 0 0];

% qpsk modulate
% mapping as defined in gr-mapper
preSyms = [];
jj = 1;
for ii=1:2:length(preSeqBits)
    bitGrp = preSeqBits(ii:ii+1);
    if(bitGrp==[0 0])
        sym = (sqrt(2)+j*sqrt(2))/2;
    elseif(bitGrp==[0 1])
        sym = (-sqrt(2)+j*sqrt(2))/2;
    elseif(bitGrp==[1 0])
        sym = (sqrt(2)-j*sqrt(2))/2;
    else
        sym = (-sqrt(2)-j*sqrt(2))/2;
    end
    preSyms(jj) = sym;
    jj = jj + 1;
end

end

