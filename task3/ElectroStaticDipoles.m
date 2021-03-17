function [Q, D] = ElectroStaticDipoles(XYZ,R, F)

n=size(R); %Количество шаров, с которыми мы работаем
n=n(2);

C=zeros(4*n); %Матрица потенциальных коэфицентов

for i=1:4*n
   for j=1:4*n
      if (i==j) && (i<=n)
         C(i, j)=1/R(i);
      elseif (i~=j) && (i<=n) && (j<=n)
          l=((XYZ(1, i)-XYZ(1, j)).^2+(XYZ(2, i)-XYZ(2, j)).^2+(XYZ(3, i)-XYZ(3, j)).^2).^0.5;
          if l>R(i)+R(j)
            C(i, j)=1/l;
          else
            error('Balls conflict detected!');
          end
      elseif (i<=n || j<=n)
          if mod(abs(i-j), n)==0
              C(i, j)=0;
          else
              l=((XYZ(1, mod(i-1, n)+1)-XYZ(1, mod(j-1, n)+1)).^2+(XYZ(2, mod(i-1, n)+1)-XYZ(2, mod(j-1, n)+1)).^2+(XYZ(3, mod(i-1, n)+1)-XYZ(3, mod(j-1, n)+1)).^2).^0.5;
              ochenyvazhnyanumber=round((max(i, j)-n-1)/n+0.5);
              g=(XYZ(ochenyvazhnyanumber, mod(i-1, n)+1)-XYZ(ochenyvazhnyanumber, mod(j-1, n)+1));
              C(i, j)=g/l^3;
          end
       else
           if ((abs(i-j)==n) || (abs(i-j)==2*n) || (abs(i-j)==3*n) || (abs(i-j)==0))
             ochenyvazhnyanumber_i=round((i-n-1)/n+0.5);
             ochenyvazhnyanumber_j=round((j-n-1)/n+0.5);
             if ochenyvazhnyanumber_i==ochenyvazhnyanumber_j
                C(i, j)=-1/(R(mod(i-1, n)+1)^3);
             else
                C(i, j)=0;
             end
           else
             l=((XYZ(1, mod(i-1, n)+1)-XYZ(1, mod(j-1, n)+1)).^2+(XYZ(2, mod(i-1, n)+1)-XYZ(2, mod(j-1, n)+1)).^2+(XYZ(3, mod(i-1, n)+1)-XYZ(3, mod(j-1, n)+1)).^2).^0.5;
             ochenyvazhnyanumber_i=round((i-n-1)/n+0.5);
             ochenyvazhnyanumber_j=round((j-n-1)/n+0.5);
             g_i=(XYZ(ochenyvazhnyanumber_i, mod(i-1, n)+1)-XYZ(ochenyvazhnyanumber_i, mod(j-1, n)+1));
             g_j=(XYZ(ochenyvazhnyanumber_j, mod(i-1, n)+1)-XYZ(ochenyvazhnyanumber_j, mod(j-1, n)+1));
             if ochenyvazhnyanumber_i==ochenyvazhnyanumber_j
                 C(i, j)=3*g_i*g_j/(l^5)-1/(l^3);
             else 
                 C(i, j)=3*g_i*g_j/(l^5);
             end
           end
      end
   end
end

D=zeros(3*n, 1);
F=[F;D];

Q=C\F;
D=[Q(n+1:2*n), Q(2*n+1:3*n), Q(3*n+1:4*n)];
Q=Q(1:n);
end

