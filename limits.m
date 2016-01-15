
function [Beval error errorF] = limits(R,i)

%%R is a time series with the columns corresponding to the species and 
%   lines to different time points. If the first column of your time series
%   labels time, remove it before applying limits (ie, do R=R(:,2:end)). 
%%i is the column of the interaction matrix we try to reconstruct.

%%Beval : the output of the function is the estimation of the "ith" line 
%         of the interaction matrix 
% error : mean of the errors made on evalution of "y" using the estimated B matrix normalized by the variance of y ("one-time-step evalutaion"). 
% errorF : idem as error but without the normalization by the variance.
listnumbkeysp=[];%list of number species kept. 

% choices to be made:
thresh = .5;% orignially put to 5 (diminish to increase precision)
r = 100; %% the number of iterations used for Bagging 


% manipulating R to put the data in the appropriated form to apply limits:
R(R == 0) = realmin;%we will have to take the log, so let's remove zero's.
sd=size(R);
N=sd(2);%number of species
Ntp=sd(1)-1;%number of time points
data= R-repmat(median(R),sd(1),1);%first N column of data matrix needed for limits
data=data(1:(sd(1)-1),:);
data(:,(sd(2)+1))=log(R(2:sd(1),i))-log(R(1:(sd(1)-1),i));


clearvars res excluded included data1 data2 test trainset testset minL locL


% variable initiation
res = zeros(N,1); %% array for storing results
error=0;
errorF=0;
listspecies=1:N; %% to construct the choices of excluded/includes species



for kkk=1:r %% we do r times the same thing to get r estimations the interaction matrix B
 
%%initialize covariates to be included/excluded from regression model
 
excluded = listspecies([1:(i-1),(i+1):end]);
included = i;
 
%%randomly partition data into training and test sets
 
 trainset = randsample(1:Ntp, round(Ntp/2));
 testset = setdiff(1:Ntp,trainset);
 
 data1 = data(trainset,:);
 data2 = data(testset,:);
 
 test = included;
 [errorEt B1t] = restrictedleastsquare(data1,data2,test); %perform the estimation
 
% loop that adds covariates(=species) as long as prediction error decreases by thresh
xxx=1;
yyy=1;

while xxx == 1 && yyy ~= N
yyy=yyy+1;
%create list of performances for all regression models including an additional covariate
  test =  [included excluded(1)];
  clearvars errorlist
  errorlist(1) =  restrictedleastsquare(data1,data2,test);
 for kk=2:(length(excluded)-1)
  test = [included excluded(kk)];
  errorlist(kk)= restrictedleastsquare(data1,data2,test);
  end
  
%sort the list so that prev[[1]] has the lowest prediction error 
[minL locL]=min(errorlist);
test = [included excluded(locL)];%we choose the test set with the smallest error
[errorEtemp B1temp errorFt]  =  restrictedleastsquare(data1, data2,test); %we re-run limits for that test set since we didn't store results
% redefine included excluded covariates to account for including new covariate 
included = test;
excluded = setxor(listspecies,test);
% if prev[[1]] improves over the current best model by thresh,  include it, otherwise exit loop 
  
  if 100.0*(errorEtemp - errorEt)/errorEt < -thresh %we keep adding species
      clearvars B1t
      errorEt=errorEtemp;% we update the error to be compared with. 
      B1t= B1temp;
      testt= test ;
      errorF=errorFt;
      
  else % we stop adding species
      xxx=0;
      listnumbkeysp=[listnumbkeysp,yyy];
  end
end
% store final regression coefficients in res 
 clearvars B
 B = zeros(N,1);
 B(test)=B1temp;

 res=[res, B];
 error=[error, errorEt];
 errorF=[errorF, errorFt];

    
 

end


% Bagging step: output median of res 
res=res(:,2:end);
Beval=median(res,2)';
error=error(2:end);
errorF=errorF(2:end);
error=median(error);
errorF=median(errorF);
listnumbkeysp;
end
%save('evaluation_b.txt','Beval','-ASCII')