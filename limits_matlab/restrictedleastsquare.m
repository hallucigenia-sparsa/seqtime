function [errorE B1 errorF]= restrictedleastsquare(inbag,outbag,test)
   % This subroutine performs the linear regression and computes the error on the test set
   %"inbag" is a matrix containing the training set. One row of the matrix looks like {x_1(t) - u_1, ..., x_N(t) - 
   % u_N, ln x_i(t+1) - ln x_i(t) }, 
    %where u_i is the median abundance of species i. ;
    %"outbag" is a matrix containing the test set. 
    %It is in the same form as inbag.;
    %"test" is a vector of integers specifying which covariates are included in the regression. I.e. if species 1,2,3 are included then test = {1,2,3}.

   %perform linear regression on training set
    
    X1 = inbag(:,test);
    y1 = inbag(:, end);
    B1 = pinv(X1)*y1;
    
    % calculate prediction error on test set 
    
    X2 = outbag(:,test);
    y2 = outbag(:,end);
    errorE = mean((y2-X2*B1).*(y2-X2*B1))/var(y2);
    errorF = mean((y2-X2*B1).*(y2-X2*B1));
    
end
  