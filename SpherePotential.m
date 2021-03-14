function [F,X,Y,P] = SpherePotential(XYZ,Q,R,r0,a,b,Dx,Dy,Nxy)
P=[0 0;0 0;0 0];
c=cross(b, a);%Вектор вдоль оси OY
if norm(c)<10^(-10)
   error('a and b should not be parallel'); 
end
P(:, 1)=a;
P(:, 2)=c;%Сформировали матрицу перехода
F=zeros(Nxy(1), Nxy(2));
Y=zeros(Nxy(1), Nxy(2));
X=zeros(Nxy(1), Nxy(2));

n=size(Q);
n=n(1);%количество шаров
for i=1:Nxy(1)
   for j=1:Nxy(2)
      X(i, j)=Dx(1)+(j-1)/(Nxy(2)-1)*(Dx(2)-Dx(1));
      Y(i, j)=Dy(1)+(i-1)/(Nxy(1)-1)*(Dy(2)-Dy(1));
      r_original=r0+P*[X(i, j);Y(i, j)];%Вектор в координатах исходной СК
      for k = 1:n
          l=norm(XYZ(:, k)-r_original);
          if l>R(k)
             F(i, j)=F(i, j)+Q(k)/l;
          else
             F(i, j)=F(i, j)+Q(k)/R(k);
          end
      end
   end
end
end

