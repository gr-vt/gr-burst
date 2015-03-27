function [outVec] = prepareCmplxVecForWrite(inVec)
% expect complex data in row or column vector form

% convert to row vectors
inVec = inVec(:);

%split it
inputDataSplit = [real(inVec) imag(inVec)];
outVec = reshape(inputDataSplit', 1, size(inputDataSplit,1)*size(inputDataSplit,2));

end
