function [Fr,Dr] = VibraStates(xyz,q,m,k) 
n=size(q);
n=n(2);
F_eq=zeros(1, 3*n); %Значение сил, действующих на i атом
C=zeros(3*n, 3*n); %Матрица dF_i/dr_j, условие на собственную моду dF/dr *u=-w^2 *u
K_l=zeros(n, n);%Равновесные значения длин связей
K_f=zeros(n, n);%Коэфиценты упрогости
M=size(k);
M=M(2);

K_el=14.39964485; %Коэфицент в законе кулона при еденицах А,еВ
Freq_=9.787151569*10^13;%Коэфицент пересчета частоты в рад/с

for i=1:M
    K_f(k(1, i), k(2, i))=k(3, i);
    K_f(k(2, i), k(1, i))=k(3, i);
    K_l(k(1, i), k(2, i))=k(4, i);
    K_l(k(2, i), k(1, i))=k(4, i);
end
for i=1:n
    for j=1:n
        if i~=j
            l=sqrt((xyz(1, i)-xyz(1, j))^2+(xyz(2, i)-xyz(2, j))^2+(xyz(3, i)-xyz(3, j))^2);
            F_eq(i)=F_eq(i)+K_el*q(i)*q(j)*(xyz(1, i)-xyz(1, j))/l^3;
            F_eq(i)=F_eq(i)-K_f(i, j)*(xyz(1, i)-xyz(1, j))*(l-K_l(i, j))/l;
            F_eq(i+n)=F_eq(i+n)+K_el*q(i)*q(j)*(xyz(2, i)-xyz(2, j))/l^3;
            F_eq(i+n)=F_eq(i+n)-K_f(i, j)*(xyz(2, i)-xyz(2, j))*(l-K_l(i, j))/l;
            F_eq(i+2*n)=F_eq(i+2*n)+K_el*q(i)*q(j)*(xyz(3, i)-xyz(3, j))/l^3;
            F_eq(i+2*n)=F_eq(i+2*n)-K_f(i, j)*(xyz(3, i)-xyz(3, j))*(l-K_l(i, j))/l;
        end
    end
    if abs(F_eq(i))>0.00001
       error('Started config. is not equilibrium');
    end
end
%Сейчас будет ужасный цикл, определяющий dF/dr для каждого атома
dr=0.000001; %dr - малое отклонение, чтобы определить dF/dr, мне лень брать производную
for i=1:n
    for j=1:n
        if i~=j
            l=sqrt((xyz(1, i)-xyz(1, j))^2+(xyz(2, i)-xyz(2, j))^2+(xyz(3, i)-xyz(3, j))^2);
            ldx=sqrt((xyz(1, i)-xyz(1, j)-dr)^2+(xyz(2, i)-xyz(2, j))^2+(xyz(3, i)-xyz(3, j))^2);
            ldy=sqrt((xyz(1, i)-xyz(1, j))^2+(xyz(2, i)-xyz(2, j)-dr)^2+(xyz(3, i)-xyz(3, j))^2);
            ldz=sqrt((xyz(1, i)-xyz(1, j))^2+(xyz(2, i)-xyz(2, j))^2+(xyz(3, i)-xyz(3, j)-dr)^2);
            
            C(i, j)=C(i, j)-K_el*q(i)*q(j)*(xyz(1, i)-xyz(1, j))/l^3+K_el*q(i)*q(j)*(xyz(1, i)-xyz(1, j)-dr)/ldx^3;
            C(i, j)=C(i, j)+K_f(i, j)*(xyz(1, i)-xyz(1, j))*(l-K_l(i, j))/l-K_f(i, j)*(xyz(1, i)-xyz(1, j)-dr)*(ldx-K_l(i, j))/ldx;
            
            C(n+i, j)=C(n+i, j)-K_el*q(i)*q(j)*(xyz(2, i)-xyz(2, j))/l^3+K_el*q(i)*q(j)*(xyz(2, i)-xyz(2, j))/ldx^3;
            C(n+i, j)=C(n+i, j)+K_f(i, j)*(xyz(2, i)-xyz(2, j))*(l-K_l(i, j))/l-K_f(i, j)*(xyz(2, i)-xyz(2, j))*(ldx-K_l(i, j))/ldx;
            
            C(2*n+i, j)=C(2*n+i, j)-K_el*q(i)*q(j)*(xyz(3, i)-xyz(3, j))/l^3+K_el*q(i)*q(j)*(xyz(3, i)-xyz(3, j))/ldx^3;
            C(2*n+i, j)=C(2*n+i, j)+K_f(i, j)*(xyz(3, i)-xyz(3, j))*(l-K_l(i, j))/l-K_f(i, j)*(xyz(3, i)-xyz(3, j))*(ldx-K_l(i, j))/ldx;
            
            
            C(i, j+n)=C(i, j+n)-K_el*q(i)*q(j)*(xyz(1, i)-xyz(1, j))/l^3+K_el*q(i)*q(j)*(xyz(1, i)-xyz(1, j))/ldy^3;
            C(i, j+n)=C(i, j+n)+K_f(i, j)*(xyz(1, i)-xyz(1, j))*(l-K_l(i, j))/l-K_f(i, j)*(xyz(1, i)-xyz(1, j))*(ldy-K_l(i, j))/ldy;
            
            C(n+i, j+n)=C(n+i, j+n)-K_el*q(i)*q(j)*(xyz(2, i)-xyz(2, j))/l^3+K_el*q(i)*q(j)*(xyz(2, i)-xyz(2, j)-dr)/ldy^3;
            C(n+i, j+n)=C(n+i, j+n)+K_f(i, j)*(xyz(2, i)-xyz(2, j))*(l-K_l(i, j))/l-K_f(i, j)*(xyz(2, i)-xyz(2, j)-dr)*(ldy-K_l(i, j))/ldy;
            
            C(2*n+i, j+n)=C(2*n+i, j+n)-K_el*q(i)*q(j)*(xyz(3, i)-xyz(3, j))/l^3+K_el*q(i)*q(j)*(xyz(3, i)-xyz(3, j))/ldy^3;
            C(2*n+i, j+n)=C(2*n+i, j+n)+K_f(i, j)*(xyz(3, i)-xyz(3, j))*(l-K_l(i, j))/l-K_f(i, j)*(xyz(3, i)-xyz(3, j))*(ldy-K_l(i, j))/ldy;
            
            
            C(i, j+2*n)=C(i, j+2*n)-K_el*q(i)*q(j)*(xyz(1, i)-xyz(1, j))/l^3+K_el*q(i)*q(j)*(xyz(1, i)-xyz(1, j))/ldz^3;
            C(i, j+2*n)=C(i, j+2*n)+K_f(i, j)*(xyz(1, i)-xyz(1, j))*(l-K_l(i, j))/l-K_f(i, j)*(xyz(1, i)-xyz(1, j))*(ldz-K_l(i, j))/ldz;
            
            C(n+i, j+2*n)=C(n+i, j+2*n)-K_el*q(i)*q(j)*(xyz(2, i)-xyz(2, j))/l^3+K_el*q(i)*q(j)*(xyz(2, i)-xyz(2, j))/ldz^3;
            C(n+i, j+2*n)=C(n+i, j+2*n)+K_f(i, j)*(xyz(2, i)-xyz(2, j))*(l-K_l(i, j))/l-K_f(i, j)*(xyz(2, i)-xyz(2, j))*(ldz-K_l(i, j))/ldz;
            
            C(2*n+i, j+2*n)=C(2*n+i, j+2*n)-K_el*q(i)*q(j)*(xyz(3, i)-xyz(3, j))/l^3+K_el*q(i)*q(j)*(xyz(3, i)-xyz(3, j)-dr)/ldz^3;
            C(2*n+i, j+2*n)=C(2*n+i, j+2*n)+K_f(i, j)*(xyz(3, i)-xyz(3, j))*(l-K_l(i, j))/l-K_f(i, j)*(xyz(3, i)-xyz(3, j)-dr)*(ldz-K_l(i, j))/ldz;
        end
    end
end
for i=1:n
   for j=1:n
      if i~=j
         for x1=0:2
             for y1=0:2
                C(i+x1*n, i+y1*n)=C(i+x1*n, i+y1*n)-C(i+x1*n, j+y1*n); 
             end
         end
      end
   end
end
C=C./dr;
C;
S=eye(3*n);
for i=1:n
   S(i, i)=1/sqrt(m(i));
   S(n+i, n+i)=1/sqrt(m(i));
   S(2*n+i, 2*n+i)=1/sqrt(m(i));
end
C=S*C*S; % Для диагонализации det(C-w^2*M)=0=>S*m*S=E, det(S*C*S-w^2*E)=0
C;

[Vectors, Diag] = eig(C);
Vectors;
Diag;
Fr = [];
Dr = [];
for i=1:3*n
    if abs(Diag(i, i))>0.001
       Fr=[Fr, Freq_*sqrt(-Diag(i, i))];
       Dr=[Dr, Vectors(:, i)];
    end
end

end

