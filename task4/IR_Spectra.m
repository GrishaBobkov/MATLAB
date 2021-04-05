function [I] = IR_Spectra(Fr,Dr,q, m, T)
%как известно d''=w^2 q u, а амплитуда Sum(m u^2 w^2/2)=E(w)=hw(1/2+1/(e^hw/kT + 1))

h=6.582119281*10^(-16);%Постоянная планка с чертой в эВ*с
k=8.617332385*10^(-5);%Постоянная больцмана
Energy_coef=1.043968445*10^(-28);%Коэффицент пересчета энергии из а.е.м*А^2*с^2 в эВ
L=size(T);
L=L(1);
n=size(q);
n=n(2);
n_freq=size(Fr);
n_freq=n_freq(2);
I=zeros(L, n_freq);

% DDr=zeros(n, n_freq);
% for i=1:n
%     for j=1:n_freq
%         DDr(i, j)=sqrt(Dr(i, j)^2+Dr(i+n, j)^2+Dr(i+2*n, j)^2);
%     end
% end
for i=1:L
    for j=1:n_freq
        %Вычисляем амплитуду колебаний.
        E=0;
        for x=1:n
            E=E+Energy_coef*m(x)*(Dr(x, j)^2+Dr(x+n, j)^2+Dr(x+2*n, j)^2)*(Fr(j)^2)/2;
        end
        ehwkt=exp(-h*Fr(j)/(k*T(i)));
        
        a=h*Fr(j)*(1/2+ehwkt/(1-ehwkt))/E;
        d=[0 0 0]; 
        for x=1:n
            d(1)=d(1)+4.8*10^(-18)*a*Dr(x, j)*q(x)*Fr(j)^2;         % В СГС, пересчет 1e*1A*c^-2=4.8 10^-18 СГСЭ
            d(2)=d(2)+4.8*10^(-18)*a*Dr(x+n, j)*q(x)*Fr(j)^2;
            d(3)=d(3)+4.8*10^(-18)*a*Dr(x+2*n, j)*q(x)*Fr(j)^2;
        end
        I(i, j)=I(i, j)+(2/(81*10^(30)))*(d(1)^2+d(2)^2+d(3)^2); 
    end
end
end

