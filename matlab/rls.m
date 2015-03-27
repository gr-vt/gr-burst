function [W,E] = rls(x,d,nord,lambda)
%RLS	Recursive Least Squares.
%--- 
%USAGE	[W,E] = rls(x,d,nord,lambda)
%
%           x    : input data to the adaptive filter.
%           d    : desired output
%           nord : number of filter coefficients
%           lambda : exponential forgetting factor
%
%     The output matrix W contains filter coefficients.
%        - The n'th row contains the filter coefficients at time n
%        - The m'th column contains the m'th filter coeff vs. time.
%        - The output vector E contains the error sequence versus time.
%
%  see also LMS and NLMS
%
%---------------------------------------------------------------
% copyright 1996, by M.H. Hayes.  For use with the book 
% "Statistical Digital Signal Processing and Modeling"
% (John Wiley & Sons, 1996).
%---------------------------------------------------------------

delta=0.001;
X=convm(x,nord);
[M,N] = size(X);
if nargin < 4,   lambda = 1.0;   end
P=eye(N)/delta;
W(1,:)=zeros(1,N);
for k=2:M-nord+1;
    z=P*X(k,:)';
    g=z/(lambda+X(k,:)*z);
    alpha=d(k)-X(k,:)*W(k-1,:).';
    W(k,:)=W(k-1,:)+alpha*g.';
    P=(P-g*z.')/lambda;
end;






