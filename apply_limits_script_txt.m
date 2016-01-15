%This scrpit computes the mean (and median) errors made on time series
%estimations for fixed N, tend, sigma.

close all 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TO FILL BEFORE RUNNING CODE                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
folderSOC='/Users/sdebuyl/Dropbox/bio/dynamics_of_microbial_communities/matlab_files/SOC_Simulations';
folderHubbell='/Users/sdebuyl/Dropbox/bio/dynamics_of_microbial_communities/matlab_files/Hubbell_Simulations';
if ~isdir('SOC_Simulations')%create a folder to save data if it doesn't exist
    mkdir('SOC_Simulations');
end
if ~isdir('Hubbell_Simulations')%create a folder to save data if it doesn't exist
    mkdir('Hubbell_Simulations');
end
%FYI : Ricker data is in folder (where your .m files are/Ricker_Noise)
maxN=10; %maximal number of species kept (we keep the most abundant ones)
maxT=500; %maximal number of time steps kept
%choose type of time series you want to analyse
fromhubbell=false; 
fromsoc=false;
generatedfromricker=false;
datafromrickermultispeciesnoise=true;
fromklemm=true;
%mix up the time series? (to test the code)
randomize=false;
%if you want to plot a reduced inferred interaciton matrix with better resolution 
plotAred=false; 
uniformd=false;%if uniformly distributed interactions

  if fromklemm==true 
      
    
    sigma=0.05;%noise in the time series used to infer the interaction matrix
    c=0.05;
    percneg=90;
    N=100;
    runnumber=1;
    fileName=strcat('klemm_notrans_',num2str(N),'_c_',num2str(c),'_clique_',num2str(5),'_percneg_',num2str(percneg),'_no_para_true',num2str(runnumber));
         
    namefiledata=strcat(fileName,'_data.txt');
    if randomize==true
            fileName=strcat(fileName,'_rand');
    end
  else
      namefiledata='test';
        if fromhubbell==true
            namefiledata='hubbell_data';
        end
        if fromsoc==true
            namefile='soc_data';
        end
  end
    
    
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END of "TO FILL BEFORE RUNNING CODE"       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 





 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% import data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

if datafromrickermultispeciesnoise==true | generatedfromricker==true
    % import our ts 
    originalDirectory=pwd;
    folder ='Ricker_Noise';   
    A=importdata('matrix.txt');
    if  wkeep(size(A),1)~=N
    return
    end
    K=importdata('carcap.txt');
    Korign=importdata('carcap.txt');
    data=importdata('data.txt');
    data=data(:,2:end);
    sd=size(data);
    tend=sd(1);
    N=sd(2);
    
  
             
  
end

if datafromrickermultispeciesnoise==false && generatedfromricker==false

    
    if fromsoc==true
%%FRANZISKA DATA 
    folder=folderSOC;
    fileName='SOI_1_adjacency_N100_c0.05_clique5_negedgeperc10_noparasites';
    
            
    originalDirectory=pwd;
    cd(folder);
    ddd=importdata('SOI_1_adjacency_N100_c0.05_clique5_negedgeperc10_noparasites_1.csv');
    cd(originalDirectory);
    data=ddd.data;
    %data=data(1:100,1:10);
    sd=size(data);
    tend=sd(1);
    N=sd(2);
    
    namefiledata='SOC_data.txt';
     
    cd(folder);         
    fID=fopen(namefiledata,'w'); % file ID ('w' means "open or create new file for Writing")
    myformat='%s\t%g\n'; 
    fprintf(fID,myformat,fileName);
    fclose(fID); 
    cd(originalDirectory);
    end

if fromhubbell==true
%           %%%LEO DATA 
    folder=folderHubbell;
    fileName='untbmigration0c02';
    
         
    originalDirectory=pwd;
    cd(folder);
    ddd=importdata(strcat(fileName,'.csv'));
    cd(originalDirectory);
    data=ddd; %(1:5000/500:5000,:);
    namefiledata='Hubbell_data.txt';  
    cd(folder);
    fID=fopen(namefiledata,'w'); % file ID ('w' means "open or create new file for Writing")
    myformat='%s\t%g\n'; 
    fprintf(fID,myformat,fileName);
    fclose(fID); 
    cd(originalDirectory);
    sd=size(data);
    tend=sd(1);
    N=sd(2);
end
%     

end


         if maxN<N
             fileName=strcat(fileName,'_',num2str(maxN));
         end
         
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% name of file to keep data about the original time series & rename
% fileName if randomize is true
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
         
         
 namefiledata=strcat(fileName,'_data.txt');
    
 if randomize==true
            fileName=strcat(fileName,'_rand');
 end         
         

         
         
         
         
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% manipulate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

if randomize==true
    data=data(randsample(1:sd(1),sd(1)),:);
end
    
 if sd(1)>maxT
        data=data(1:round(sd(1)/200):sd(1),:);
    end
    
    K=mean(data);
    data=data(:,K>0);
    K=mean(data);
    %reduce number of species if necessary
    [i a]=sort(sort(K));
    [i b]=sort(K);
    mapping=b(a);
    sd=size(data);
    
    
   
    
    %data=data(1:250:5000,mapping((sd(2)-39):sd(2)));
    if sd(2)>maxN
        %data=data(1:20:5000,:);
        data=data(:,mapping((sd(2)-maxN+1):sd(2)));
    sd=size(data);
    if generatedfromricker==true | datafromrickermultispeciesnoise==true;
        A=A(mapping((sd(2)-maxN+1):sd(2)),mapping((sd(2)-maxN+1):sd(2)));
    end
    end
    
    sd=size(data);
    tend=sd(1);
    N=sd(2);
    Z=sum(data(tend,:));
    data=data/Z;
    
   
    K=mean(data);

    
    
    
 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN LIMITS                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

[Best error errorf] = limits(data,1);


for i=2:sd(2)
    i
[Beval errort errorft] = limits(data,i);
Best=[Best;Beval];
error=[error errort];
errorf=[errorf errorft];
end

%realBest=Best/Z;

%Best(logical(eye(size(Best)))) = -1;% impose that it is -1 on the diagonal. 

%check linear stability of infered matrix & remove add minus sign to
%positive real parts of eigenvalues. 
eigenvalBest=eig(Best); 
if max(real(eigenvalBest))>=0
    fprintf('linear unstable interaction matrix')
    
    [U,T]=schur(Best,'complex');

    rediag=real(diag(T));
    imdiag=imag(diag(T));
    indicesP=rediag>0;

    listind=1:N;

    for k=listind(indicesP)
    T(k,k)=-real(T(k,k))+i*imag(T(k,k));
    end

    Best=U*T*U';
    max(max(imag(Best)))
    Best=real(Best);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE ERRORS                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%compute time series with the same initial condition and same interaction
%matrix and compute the "difference" with the original time series.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

x=data(1,:)';
y=data(1,:)';
Restim=x';
RestimLoc=y';

  RestimLoc2=[data(1,:); data(2,:)];
  RestimLoc3=[data(1,:); data(2,:); data(3,:)];
  RestimLoc4=[data(1,:); data(2,:); data(3,:);data(4,:)];
  RestimLoc5=[data(1,:); data(2,:);data(3,:);data(4,:);data(5,:)];

Kestim=mean(data(round(tend*3/4):end,:))';
Kestim=mean(data)';
BestZ=Best;

for t=2:tend
    x=x.*exp(BestZ*(x-Kestim)); %lognrnd(0,0.05,[N,1]).*
    Restim=[Restim; x'];
    tm1=t-1;
    
    estimyt=data(tm1,:)'.*exp(BestZ*(data(tm1,:)'-Kestim));%estimation at time t (based on t-1 of original ts at time t). 
    RestimLoc=[RestimLoc; estimyt'];
   
    estimyt2=estimyt.*exp(BestZ*(estimyt-Kestim));%estimation at time t+1 base on orginal ts at time t. 
    RestimLoc2=[RestimLoc2; estimyt2'];
    
    estimyt3=estimyt2.*exp(BestZ*(estimyt2-Kestim));
        RestimLoc3=[RestimLoc3; estimyt3'];
    estimyt4=estimyt3.*exp(BestZ*(estimyt3-Kestim));
         RestimLoc4=[RestimLoc4; estimyt4'];
    estimyt5=estimyt4.*exp(BestZ*(estimyt4-Kestim));
        RestimLoc5=[RestimLoc5; estimyt5'];
     
end




KestimR=repmat(Kestim,1,sd(1)-1);
% 
% errortsloc=0;
% errortsloc2=0;
% errortsloc3=0;
% errortsloc4=0;
% errortsloc5=0;
% 
% for kk=1:N
%     temp=corrcoef(data(2:tend,kk),RestimLoc(2:tend,kk));
%     errortsloc=errortsloc+temp(1,2);
%     temp=corrcoef(data(3:tend,kk),RestimLoc2(3:tend,kk));
%     errortsloc2=errortsloc2+temp(1,2);
%      temp=corrcoef(data(4:tend,kk),RestimLoc3(4:tend,kk));
%     errortsloc3=errortsloc3+temp(1,2);
%      temp=corrcoef(data(5:tend,kk),RestimLoc4(5:tend,kk));
%     errortsloc4=errortsloc4+temp(1,2);
%     temp=corrcoef(data(6:tend,kk),RestimLoc5(6:tend,kk));
%     errortsloc5=errortsloc5+temp(1,2);
% end
% errortsloc=errortsloc/N;
% errortsloc2=errortsloc2/N;
% errortsloc3=errortsloc3/N;
% errortsloc4=errortsloc4/N;
% errortsloc5=errortsloc5/N;

%to be compared with autocorrelation
% 
% autocorr=0;
% autocorr1=0;
% autocorr2=0;
% autocorr3=0;
% autocorr4=0;
% for kk=1:N
%      temp=corrcoef(data(1:tend,kk),data(1:tend,kk));
%     autocorr=autocorr+temp(1,2);
%      temp=corrcoef(data(1:tend-1,kk),data(2:tend,kk));
%     autocorr1=autocorr1+temp(1,2);
%     temp=corrcoef(data(1:tend-2,kk),data(3:tend,kk));
%     autocorr2=autocorr2+temp(1,2);
%      temp=corrcoef(data(1:tend-3,kk),data(4:tend,kk));
%     autocorr3=autocorr3+temp(1,2);
%     temp=corrcoef(data(1:tend-4,kk),data(5:tend,kk));
%     autocorr4=autocorr4+temp(1,2);
% end
% 
% autocorr=autocorr/N;
% autocorr1=autocorr1/N;
% autocorr2=autocorr2/N;
% autocorr3=autocorr3/N;
% autocorr4=autocorr4/N;


%initialize lists of errors with reduced number of species



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% error on interaction matrix (only if data from ricker)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

if generatedfromricker==true | datafromrickermultispeciesnoise==true
    tempA=removediag(A);% we remove the diagonal before computing correlation. 
    tempB=removediag(Best/Z);
    temperroronA=corrcoef(tempA,tempB);
    erroronA=temperroronA(1,2);
    temp=corrcoef(A,Best);
    corcoefABred=temp(1,2);
    errorAred=erroronA;%initialize erroronA list 
    lserrorABred=sum(sum((A-Best/Z).*(A-Best/Z)))/N^2;
end

sizered=N;

listindices=1:N;

% errortslocred=errortsloc;
% errortsloc2red=errortsloc2;
% errortsloc3red=errortsloc3;
% errortsloc4red=errortsloc4;
% 
% autocorr1red=autocorrts(data,listindices,1);
% autocorr2red=autocorrts(data,listindices,2);
% autocorr3red=autocorrts(data,listindices,3);
% autocorr4red=autocorrts(data,listindices,4);

errortslocred=0;
errortsloc2red=0;
errortsloc3red=0;
errortsloc4red=0;

autocorr1red=0;
autocorr2red=0;
autocorr3red=0;
autocorr4red=0;

percOKlist=[.001 .01 .1 .15 .2 .25 .3 .35 .4 .45 .5 .55 .6 .65 .7 .75 .8 .85 .9 .95 100];

            for percOK=percOKlist

                indicesOK=mean(data)> percOK*max(mean(data));
                Bred=BestZ(indicesOK,indicesOK);
                
                if generatedfromricker==true | datafromrickermultispeciesnoise==true
                    Ared=A(indicesOK,indicesOK);
                end
                
                sB=size(Bred);
                sizeredT=sB(1);
                
                
                if (generatedfromricker==true | datafromrickermultispeciesnoise==true) && sizeredT>1
                    sizered=[sizered sizeredT];
                    tempA=removediag(Ared);% we remove the diagonal before computing correlation. 
                    tempB=removediag(Bred/Z);
                    temperroronA=corrcoef(tempA,tempB);
                    erroronA=temperroronA(1,2);
                    temp=corrcoef(Ared,Bred/Z);
                    corcoefAB=temp(1,2);
                    lserrorAB=sum(sum((Ared-Bred/Z).*(Ared-Bred/Z)))/sizeredT^2;
                    
                errorAred=[errorAred erroronA]; 
                corcoefABred=[corcoefABred corcoefAB];
                lserrorABred=[lserrorABred lserrorAB];
                else
                    sizered=[sizered sizeredT];
                end
                
                if sizeredT>1
                    errortsloc=0;
                    errortsloc2=0;
                    errortsloc3=0;
                    errortsloc4=0;
                    for kk=listindices(indicesOK)
                        temp=corrcoef(data(1:tend,kk),RestimLoc(1:tend,kk));
                        errortsloc=errortsloc+temp(1,2);
                         temp=corrcoef(data(2:tend,kk),RestimLoc2(2:tend,kk));
                        errortsloc2=errortsloc2+temp(1,2);
                          temp=corrcoef(data(3:tend,kk),RestimLoc3(3:tend,kk));
                        errortsloc3=errortsloc3+temp(1,2);
                        temp=corrcoef(data(4:tend,kk),RestimLoc4(4:tend,kk));
                        errortsloc4=errortsloc4+temp(1,2);
                    end
                    
                    errortsloc=errortsloc/sizeredT;
                    errortslocred=[errortslocred, errortsloc];
                   
                    errortsloc2=errortsloc2/sizeredT;
                    errortsloc2red=[errortsloc2red, errortsloc2];
                    
                    errortsloc3=errortsloc3/sizeredT;
                    errortsloc3red=[errortsloc3red, errortsloc3];
                    
                    errortsloc4=errortsloc4/sizeredT;
                    errortsloc4red=[errortsloc4red, errortsloc4];
                    
                    autocorr1=autocorrts(data,listindices(indicesOK),1);
                    autocorr1red=[autocorr1red,autocorr1];
                    
                    autocorr2=autocorrts(data,listindices(indicesOK),2);
                    autocorr2red=[autocorr2red,autocorr2];
                    
                    autocorr3=autocorrts(data,listindices(indicesOK),3);
                    autocorr3red=[autocorr3red,autocorr3];
                    
                    autocorr4=autocorrts(data,listindices(indicesOK),4);
                    autocorr4red=[autocorr4red,autocorr4];
                    
                end

% % %                 xred=data(1,indicesOK)';
% % %                 Restimred=xred';
% % %                 Kestimred=Kestim(indicesOK);
% % %                 for t=2:tend
% % %                     xred=lognrnd(0,sigma,[sB(1),1]).*xred.*exp(Bred*(xred-Kestimred)); %lognrnd(0,0.05,[N,1]).*
% % %                     Restimred=[Restimred; xred'];
% % %                 end
% % %                 
% % %                 errortsredT=mean(diag(corr(data(:,indicesOK),Restimred)));
% % %                 
% % %                 
% % %                 errortsred=[errortsred errortsredT];

            end
            
            
errortslocred=errortslocred(2:end);
errortsloc2red=errortsloc2red(2:end);
errortsloc3red=errortsloc3red(2:end);
errortsloc4red=errortsloc4red(2:end);

autocorr1red=autocorr1red(2:end);
autocorr2red=autocorr2red(2:end);
autocorr3red=autocorr3red(2:end);
autocorr4red=autocorr4red(2:end);

sizered=sizered(2:end);
            
sizeT=size(sizered(sizered>1));
sizeT=sizeT(2);





%%IF YOU WANT TO HAVE A LOOK at exploding cases, add this line :
%listnoexplosion= setdiff(1:nbrruns,listnoexplosion);














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTS (PNG)                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% figures correlation coef between A and its estimation (as a function of number of species kept) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            if generatedfromricker==true | datafromrickermultispeciesnoise==true
            h=figure;
            clf
            plot(sizered(sizered>1),errorAred(sizered>1),'b.--')
            hold on
            plot(sizered(sizered>1),corcoefABred(sizered>1),'r.--')
            hold on
            plot(sizered(sizered>1),lserrorABred(sizered>1),'g.--')
            legend('corrcoef between A and its estimation (without diagonal)','corrcoef between A and its estimation ','mean squared error made on each element of the interaction matrix ','Location','southwest');

            title('Comparison of the interaction matrix and its estimation','fontsize',16)
            axis([1 N -1 1])
            xlabel('Number of species included to compute correlation','fontsize',12); 
            ylabel('R coeff','fontsize',12);
            cd(folder);
            figNamej=strcat(fileName,'_errorA.png');
            saveas(h,figNamej,'png');
            cd(originalDirectory);
            hold off
            end
            
            h=figure;
            clf
            plot(sizered(sizered>1),errortslocred(sizered>1),'b.--')
            hold on
             plot(sizered(sizered>1),errortsloc2red(sizered>1),'r.--')
             hold on
               plot(sizered(sizered>1),errortsloc3red(sizered>1),'k.--')
             hold on
               plot(sizered(sizered>1),errortsloc4red(sizered>1),'m.--')
             hold on
             plot(sizered(sizered>1),autocorr1red(sizered>1),'b*')
             hold on
             plot(sizered(sizered>1),autocorr2red(sizered>1),'r*')
             hold on
             plot(sizered(sizered>1),autocorr3red(sizered>1),'k*')
             hold on
             plot(sizered(sizered>1),autocorr4red(sizered>1),'m*')
             legend('corr ts estim ts 1 step','corr ts estim ts 2 steps','3','4','auto corr 1 step','auto corr 2 step','3','4','Location','southwest');
            title('Correlation between ts and its local estimation','fontsize',16)
            xlabel('Number of species included to compute correlation','fontsize',12); 
            ylabel('R coeff','fontsize',12);
            axis([1 N -1 1])
            cd(folder);
            figNamej=strcat(fileName,'_errorts.png');
            saveas(h,figNamej,'png');
            cd(originalDirectory);
            hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% figures of 5 species ts and the corresponding predictions 
% ONLY FOR VISUAL CHECKING OF THE CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % h=figure;
% % % clf
% % %     for nn=1:round(N/5):N
% % %     randcolor = [rand() rand() rand()];
% % %     plot(0:(tend-1),data(:,nn),'Color',randcolor);
% % %     hold on
% % %     plot(0:(tend-1),RestimLoc(:,nn),'Color',randcolor,'LineStyle', '--');%plot(Restim(:,1),[Restim(:,nn+1),R(:,nn)],'Color',randcolor);
% % %     hold on
% % %     end
% % % title('Original ts and loc estimation one time point ahead (5 species)')
% % %     xlabel('Time','fontsize',16); 
% % %     ylabel('X_i','fontsize',16);
% % %     %legend(legendString);
% % %     figNamej=strcat(fileName,'_ts_and_estim.png');
% % %     cd(folder);
% % %     saveas(h,figNamej,'png');
% % %     cd(originalDirectory);
% % %     hold off

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% figures of 5 species ts and the corresponding LOCAL  predictions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

selectsp=randperm(N);   
listtimeplot=1:round(tend/10):tend; % choose some times to plot local estimations
h=figure;
clf
    for nn=selectsp(1,1:5)
    plot(0:tend-1,data(:,nn))
    hold on
    
    for kk=listtimeplot(1:end-2)
    plot([0 1 2 3 4]+kk,[data(kk+1,nn) RestimLoc(kk+2,nn) RestimLoc2(kk+3,nn) RestimLoc3(kk+4,nn) RestimLoc4(kk+5,nn)],'*r-')
    hold on
    end
    end
    xlabel('Time','fontsize',16); 
    ylabel('X_i','fontsize',16);
    title('ts and local estimation 5 time points ahead');
    
    %legend(legendString);
    figNamej=strcat(fileName,'_ts_and_estimLOC.png');
    cd(folder);
    saveas(h,figNamej,'png');
    cd(originalDirectory);
    hold off
    
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% figures of TS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=figure;
clf
    for nn=1:N
    randcolor = [rand() rand() rand()];
    plot(0:(tend-1),data(:,nn),'Color',randcolor);
    hold on
    end

    xlabel('Time','fontsize',16); 
    ylabel('X_i','fontsize',16);
    title('original ts');

    %legend(legendString);
    figNamej=strcat(fileName,'_ts.png');
    cd(folder);
    saveas(h,figNamej,'png');
    cd(originalDirectory);
    hold off
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% figures of estimated TS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=figure;
clf
    for nn=1:N
    randcolor = [rand() rand() rand()];
    plot(0:(tend-1),Restim(:,nn),'Color',randcolor);
    hold on
    end
    title('estimated ts');
    xlabel('Time','fontsize',16); 
    ylabel('X_i','fontsize',16);
    %legend(legendString);
    figNamej=strcat(fileName,'_ts_estim.png');
    cd(folder);
    saveas(h,figNamej,'png');
    cd(originalDirectory);
    hold off


    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% figures of interaction matrix & predicted inter matrix (if applicable)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if generatedfromricker==true| datafromrickermultispeciesnoise==true

    figNamej=strcat(fileName,'_A.png');
    h=figure;
    h1=subplot(1,2,1);
    set(h1,'Position',[0.1 0.3 0.25*4/3 0.325*4/3]);
    plotmatrix(A)
    h2=subplot(1,2,2);
    set(h2,'Position',[0.5 0.3 0.25*4/3 0.325*4/3]);
    plotmatrix(Best/Z)
%note: pour la position les 2 premiers chiffres sont (x,y) de l'origine du
%plot et les deux derniers donne la taille max en x et y. 
    cd(folder);
    saveas(h,figNamej,'png');
    cd(originalDirectory);
    
% % %     
% % %     figNamej=strcat(fileName,'_Ared.png');
% % %     h=figure;
% % %     h1=subplot(1,2,1);
% % %     set(h1,'Position',[0.1 0.3 0.25*4/3 0.325*4/3]);
% % %     plotmatrix(Ared)
% % %     h2=subplot(1,2,2);
% % %     set(h2,'Position',[0.5 0.3 0.25*4/3 0.325*4/3]);
% % %     plotmatrix(Bred/Z)
% % % %note: pour la position les 2 premiers chiffres sont (x,y) de l'origine du
% % % %plot et les deux derniers donne la taille max en x et y. 
% % %     cd(folder);
% % %     saveas(h,figNamej,'png');
% % %     cd(originalDirectory);
    
    
else
    h=figure;
    plotmatrix(BestZ)
    figName=strcat(fileName,'_Aestim.png');
    cd(folder);
    saveas(h,figName,'png');
    cd(originalDirectory);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE ALL DATA IN TXT FILES                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




if generatedfromricker==true | datafromrickermultispeciesnoise==true
    AName=strcat(fileName,'_A.txt');
    carcapName=strcat(fileName,'_carcap.txt');
end
    AestimNamej=strcat(fileName,'_Aestim.txt');
    tsNamej=strcat(fileName,'_ts.txt');
    estimtsNamej=strcat(fileName,'_ts_estim.txt');
    
 if generatedfromricker==true | datafromrickermultispeciesnoise==true
     
    temp=A;
    cd(folder);
    save(AName,'temp','-ASCII');
    cd(originalDirectory);
 
    
    
    cd(folder);
    save(carcapName,'temp','-ASCII');
    cd(originalDirectory);
    
 end
    temp=Best;
    cd(folder);
    save(AestimNamej,'temp','-ASCII');
    cd(originalDirectory);
    temp=data.*Z;
    cd(folder);
    save(tsNamej,'temp','-ASCII');
    cd(originalDirectory);
    temp=Restim;
    cd(folder);
    save(estimtsNamej,'temp','-ASCII');
    cd(originalDirectory);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WRITE RESULTS IN HTML FORMAT                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fileNameHTML=strcat(fileName,'HTML','.txt');


cd(folder);
fID=fopen(fileNameHTML,'w'); % file ID ('w' means "open or create new file for Writing")
myformat='%s\t%g\n'; 

    
fprintf(fID,'<tr>');
if randomize==true
fprintf(fID,'<td valign="top"> random mixing </a></td>',fileName);
else

            if fromklemm==true
                fprintf(fID,'<td valign="top"> klemm <a href="%s"> info </a></td>',namefiledata);
            else 
                    if fromsoc==true
                    fprintf(fID,'<td valign="top"> soc <a href="%s"> info </a></td>',namefiledata);
                    else 
                        if fromhubbell==true
                             fprintf(fID,'<td valign="top"> hubbell <a href="%s"> info </a></td>',namefiledata);
                        else
                             fprintf(fID,'<td valign="top"> .. <a href="%s"> info </a></td>',namefiledata);
                        end
                    end
            end

            
end

%fprintf(fID,myformat,'</b></td>');
fprintf(fID,'\n');
fprintf(fID,myformat,'<td valign="top">');
fprintf(fID,myformat,num2str(N));
fprintf(fID,myformat,'</td>');
fprintf(fID,'\n');
fprintf(fID,myformat,'<td valign="top"> uniform from [0,1] </td>');
fprintf(fID,'\n');
if uniformd==true & (generatedfromricker==true | datafromrickermultispeciesnoise==true)
    fprintf(fID,'<td valign="top"> uniform from [%g,%g] </td>',intstr(1),intstr(2));
else
    fprintf(fID,myformat,'<td valign="top"> see info </td>');
end
fprintf(fID,'\n');
fprintf(fID,myformat,'<td valign="top">-1</td>');
fprintf(fID,'\n');
fprintf(fID,myformat,'<td valign="top">');
if exist('spars')==true
fprintf(fID,myformat,num2str(spars));
else
   fprintf(fID,myformat,'see info'); 
end
fprintf(fID,myformat,'</td>');
fprintf(fID,'\n');



if generatedfromricker==true | datafromrickermultispeciesnoise==true

fprintf(fID,myformat,'<td valign="top">');
fprintf(fID,myformat,num2str(sigma));
fprintf(fID,myformat,'</td>');
fprintf(fID,'\n');

fprintf(fID,'<td valign="top"><a href="%s_carcap.txt">carcap.txt</a> </td>',fileName);
fprintf(fID,'\n');
fprintf(fID,'<td valign="top"><a href="%s_A.txt">A.txt</a></td>',fileName);
fprintf(fID,'\n');
else
fprintf(fID,myformat,'<td valign="top">');
fprintf(fID,myformat,'.');
fprintf(fID,myformat,'</td>');
fprintf(fID,'\n');
fprintf(fID,myformat,'<td valign="top"> not appl </td>');
fprintf(fID,'\n');
fprintf(fID,myformat,'<td valign="top"> not appl </td>');
fprintf(fID,'\n');
end

fprintf(fID,'<td valign="top"><a href="%s_ts.txt">ts.txt</a>[<a href="%s_ts.png">png</a>]</td>',fileName,fileName);
fprintf(fID,'\n');
fprintf(fID,'<td valign="top"><a href="%s_Aestim.txt">Aestim.txt</a></td>',fileName);
fprintf(fID,'\n');
%fprintf(fID,'<td valign="top"><a href="%s_ts_estim.txt">tsestim.txt</a>[<a
%href="%s_ts_estim.png">png</a>]</td>',fileName,fileName);  %% useless
%picture to remove
fprintf(fID,'<td valign="top"><a href="%s_ts_estim.txt">tsestim.txt</a>[<a href="%s_ts_estim.png">png</a>]</td>',fileName,fileName);
fprintf(fID,'\n');
fprintf(fID,'<td valign="top"> [<a href="%s_ts_and_estimLOC.png">png</a>]</td>',fileName);
fprintf(fID,'\n');

if datafromrickermultispeciesnoise==true | generatedfromricker==true
fprintf(fID,'<td valign="top"><a href="%s_A.png">A_and_Aestim.png</a></td>',fileName);
fprintf(fID,'\n');
        if randomize==false
        fprintf(fID,'<td valign="top" bgcolor="#ddeeef"> %1.2g <a href="%s_errorA.png">errorA.png</a> </td>', corcoefABred(1),fileName);    
        else
        fprintf(fID,'<td valign="top"> %1.2g <a href="%s_errorA.png">errorA.png</a>  </td>',  corcoefABred(1),fileName);
        end
else

fprintf(fID,myformat,'<td valign="top"> not appl </td>');
fprintf(fID,'\n');
fprintf(fID,myformat,'<td valign="top"> not appl </td>');

end


fprintf(fID,'\n');
%fprintf(fID,'<td valign="top"> %1.2g </td>', errorts);
%fprintf(fID,'\n');
if randomize==true

fprintf(fID,'<td valign="top"> %1.2g <a href="%s_errorts.png">errorts.png</a></td>', errortslocred(2),fileName);

else
fprintf(fID,'<td valign="top" bgcolor="#fde482"> %1.2g <a href="%s_errorts.png">errorts.png</a> </td>', errortslocred(2),fileName);    
end
fprintf(fID,'</tr>');
%fprintf(fID,'<td valign="top"> %1.2g </td> </tr>', median(error));


fclose(fID); 

cd(originalDirectory);


%     old stuff
%     f = @(x) sum(sum((A-x*Best).*(A-x*Best)))/(N^2-N);
%     'check that number is close to one'
%     fminbnd(f,1,100)

if plotAred==true
    
temps=sizered>70/100*N;
tempe=errortslocred>.8;

list2=find(temps+tempe==2);


indicesOK=mean(data)> percOKlist(list2(end))*max(mean(data));
                Bred=BestZ(indicesOK,indicesOK);
Ared=A(indicesOK,indicesOK)


% % %     
    figNamej=strcat(fileName,'_Ared.png');
    h=figure;
    h1=subplot(1,2,1);
    set(h1,'Position',[0.1 0.3 0.25*4/3 0.325*4/3]);
    plotmatrix(Ared)
    h2=subplot(1,2,2);
    set(h2,'Position',[0.5 0.3 0.25*4/3 0.325*4/3]);
    plotmatrix(Bred/Z)
%note: pour la position les 2 premiers chiffres sont (x,y) de l'origine du
%plot et les deux derniers donne la taille max en x et y. 
    cd(folder);
    saveas(h,figNamej,'png');
    cd(originalDirectory);
end



% to compute least square error
% % % LSloc1=0;
% % % for kk=1:N
% % %    
% % %    LSloc1=LSloc1+sum((data(2:tend,kk)-RestimLoc(2:tend,kk)).^2);
% % % end
% % % LSloc1=LSloc1/N/tend;
% % % 
% % % 
% % % LSloc2=0;
% % % for kk=1:N
% % %     
% % %    LSloc2=LSloc2+sum((data(3:tend,kk)-RestimLoc2(3:tend,kk)).^2);
% % % end
% % % LSloc2=LSloc2/N/tend;
% % % 
% % % 
% % % LSloc3=0;
% % % for kk=1:N
% % %     
% % %    LSloc3=LSloc3+sum((data(4:tend,kk)-RestimLoc3(4:tend,kk)).^2);
% % % end
% % % LSloc3=LSloc3/N/tend;
% % % 
% % % 
% % % LSloc4=0;
% % % for kk=1:N;
% % %     
% % %    LSloc4=LSloc4+sum((data(5:tend,kk)-RestimLoc4(5:tend,kk)).^2);
% % % end
% % % LSloc4=LSloc4/N/tend;
% % % 
% % % LSloc5=0;
% % % for kk=1:N
% % %     
% % %    LSloc5=LSloc5+sum((data(6:tend,kk)-RestimLoc5(6:tend,kk)).^2);
% % % end
% % % LSloc5=LSloc5/N/tend;
% % % [LSloc1 LSloc2 LSloc3 LSloc4 LSloc5]