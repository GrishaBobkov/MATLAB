function [P,sgP] = NonLinApproximator(y, r, fun, M)
P=zeros(1, M);
N=size(y);
N=N(2);
% 
%  E=Err(y, r, fun, P)
%  E=Err(y, r, fun, P+[0 0 0.01])
%  R=df(y, r, fun, P, M)

n=0;
R=1;

 while n<10000 && norm(R)>0.0001
     R=df(y, r, fun, P, M);
     lambda=(R*R')/(R*d2f(y, r, fun, P, M)*R');
     P=P-lambda*R;
     n=n+1;
%      Err(y, r, fun, P)
 end

sgP=zeros(1, M);
Er=Err(y, r, fun, P);
H=d2f(y, r, fun, P, M);
for i=1:M
   sgP(i)=sqrt(Er/H(i, i));
end
end

function [Ee] = Err(y, r, fun, P)
Ee=0;
N=size(y);
N=N(2);
for i=1:N
   Ee=Ee+(y(i)-fun(r(:, i), P)).^2;
end
% Ee=P(1)^2+P(2)^2+(P(1)+1)*(P(2)+1);
end

function [Ee] = df(y, r, fun, P, M)
ddd=10^(-5);%Разумеется только для f \approx 1
Ee=zeros(1, M);
for i=1:M
    P2=P;
    P2(i)=P2(i)+ddd;
    Ee(i)=(Err(y, r, fun, P2)-Err(y, r, fun, P))/ddd;
end
end

function [Ee] = d2f(y, r, fun, P, M)
ddd=10^(-5);%Разумеется только для f \approx 1
Ee=zeros(M, M);
for i=1:M
    P2=P;
    P2(i)=P2(i)+ddd;
    Ee(i, :)=(df(y, r, fun, P2, M)-df(y, r, fun, P, M))/ddd;
end
end


