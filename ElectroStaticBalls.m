function [Q] = ElectroStaticBalls(XYZ,R, F)
n=size(R); %Количество шаров, с которыми мы работаем
n=n(2);

C=zeros(n); %Матрица потенциальных коэфицентов

for i=1:n
   for j=1:n
      if i==j
         C(i, j)=1/R(i);
      else
          l=((XYZ(1, i)-XYZ(1, j)).^2+(XYZ(2, i)-XYZ(2, j)).^2+(XYZ(3, i)-XYZ(3, j)).^2).^0.5;
          if l>R(i)+R(j)
            C(i, j)=1/l;
          else
            error('Balls conflict detected!');
          end
      end
   end
end

Q=C\F;

end

