function [Answer] = KronigPenney(k0, m0, a0, b0, U00, Emax0)
%KronigPenney 
global h
h=1.054571817*10.^(-27)/((1.602176634*10.^(-12)*1*10.^(-14)).^0.5);%h(безразмерная)=h/(1эВ * 1г * 1нм^2)^0.5(все в СГС)
global k
global a
global b
global m
global U0
global Emax
a=a0;
b=b0;
m=m0;
U0=U00;
Emax=Emax0;


%cos(k.*(a+b))=cos(mu.*a).*cos(l.*b)-(mu.^2-l.^2)./(2*mu.*l).*sin(mu.*a).*sin(l.*b)
%mu.^2=2*m.*E/h
%l.^2=2*m.*(E-U)/h
%

f1=@f;
a_final=[];
current_length=0;
for kk = 1:length(k0)
    k=k0(kk);
    a1=[];
    for ii = U0+0.0001:0.05:30
        a2=fzero(f1, ii);
        flag=true;
        for jj = 1:length(a1)
            if abs(a1(jj)-a2)<0.0001
                flag=false;
            end
        end 
        if (flag) && (a2<Emax0) 
            a1=[a1, a2];
        end
    end
    a1=-sort(-a1);
    a1=a1';
    if kk==1
        a_final=a1;
        current_length=length(a1);
    else
        current_length;
        if current_length<length(a1)
            a_final=[a_final; NaN.*zeros(length(a1)-current_length, kk-1)];
            current_length=length(a1);
        else
            a1=[a1; NaN.*zeros(-length(a1)+current_length, 1)];
        end
        a_final;
        a1;
        a_final=[a_final, a1];
    end
    
end

%     k=k0;
%     a1=[];
%     for ii = 0.01:0.01:(2.*Emax0)
%         a2=fzero(f1, ii);
%         flag=true;
%         for jj = 1:length(a1)
%             if abs(a1(jj)-a2)<0.0001
%                 flag=false;
%             end
%         end 
%         if (flag) && (a2<0) 
%             a1=[a1, a2];
%         end
%     end
%     a_final=a1;
         
Answer=a_final;
end

function q = f(E)
global k
global a
global b
global m
global U0
global h
mu=(2*m.*E/h.^2).^0.5;
l=(2*m.*(E-U0)/h.^2).^0.5;
    q=cos(mu.*a).*cos(l.*b)-(mu.^2-l.^2)./(2*mu.*l).*sin(mu.*a).*sin(l.*b)-cos(k.*(a+b));
end
