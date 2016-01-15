r=2;  % scale factor (if r>2 => complex dynamics (oscillations, chaos), esp. if N large)
N=100;
spars=.02;
A1=-rand(N);           % interaction matrix
S1=rand(N)<spars;        % sparsity
A1=r*A1.*S1;             % interaction matrix (sparse)
%A=zeros(N);          % no interactions (in this case, X tend to K)
A1(logical(eye(size(A1)))) = -1;   % values on the diagonals

A2=-rand(N);           % interaction matrix
S2=rand(N)<spars;        % sparsity
A2=r*A2.*S2;             % interaction matrix (sparse)
%A=zeros(N);          % no interactions (in this case, X tend to K)
A2(logical(eye(size(A2)))) = -1;   % values on the diagonals

erroronA=sum(sum((A1-A2).*(A1-A2)))/N^2