function y = autocorrts(data,vectN,shift)
%autocorrelation of time series with selected species
autocorr1=0;
sizevectN=size(vectN);
Ntemp=sizevectN(2);
for kk=vectN
    temp=corrcoef(data(1:end-shift,kk),data(1+shift:end,kk));
    autocorr1=autocorr1+temp(1,2);
    
end
y=autocorr1/Ntemp;

end

