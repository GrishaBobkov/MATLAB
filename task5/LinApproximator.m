function [P,sgP] = LinApproximator(y,r,funcs)
N=size(y);
N=N(2);
K=size(r);
K=K(1);
M=size(funcs);
M=M(1);

Y=zeros(M, 1);
G=zeros(M, M);
for i=1:M
   for j=1:N
       Y(i)=Y(i)+y(j).*funcs{i}(r(:,j));
   end
end

for i=1:M
   for j=1:M
       for k=1:N
           G(i, j)=G(i, j)+funcs{i}(r(:,k)).*funcs{j}(r(:,k));
       end
   end
end
P=G\Y;

Err=0;
for i=1:N
    d=0;
    for j=1:M
       d=d+P(j).*funcs{j}(r(:,i));
    end
    Err=Err+(y(i)-d).^2;
end
Err=Err/N;
for i=1:M
   sgP(i)=sqrt(Err/(2*G(i, i)));  
end
end

