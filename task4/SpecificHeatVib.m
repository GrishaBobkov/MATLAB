function [Cv] = SpecificHeatVib(Fr,T)
h=6.582119281*10^(-16);%Постоянная планка с чертой в эВ*с
k=8.617332385*10^(-5);%Постоянная больцмана
n=size(Fr);
n=n(2);
Cv=0.*T;
for i=1:n
    hwkt=h.*Fr(i)./(k.*T);
    
    Cv=Cv+k.*(hwkt.^2).*(exp(-hwkt))./((1-exp(-hwkt)).^2);
end

end

