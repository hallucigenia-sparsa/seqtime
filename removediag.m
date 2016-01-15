function c = removediag(A)
% remove the diagonal part of a matrix and transform it into a vector. 
%   Detailed explanation goes here


    minA=min(min(A));
    sA=size(A);
    N=sA(1);
    tempA=reshape(A+ eye(N)*(10*minA),N*N,1);
    c=tempA(tempA>2*minA);

end



