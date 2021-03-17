% xyz=[0 0 1 1; 0 1 0 1; 0 0 0 0];
% r=[0.3 0.3 0.3 0.3];
% f=[10; 10; 10; 10];
% [F, X, Y, p]=SpherePotential(xyz, ElectroStaticBalls(xyz,r, f), r, [0;0;0], [1;0;0], [0;1;0], [-0.5 1.5], [-0.5 1.5], [2000 2000]);
% figure; hold on; grid on; mesh(X,Y,F); 
% 
% xyz=[0 0 1 1; 0 1 0 1; 0 0 0 0];
% r=[0.3 0.3 0.3 0.3];
% f=[10; 10; 10; 10];
% [Q, D]=ElectroStaticDipoles(xyz, r, f);
% [F, X, Y, p]=SphereDipPotential(xyz, Q, D, r, [0;0;0], [1;0;0], [0;1;0], [-0.5 1.5], [-0.5 1.5], [2000 2000]);
% figure; hold on; grid on; mesh(X,Y,F); 


load 'Lin_13_sin_pot.mat'
XYZ=XYZ';
R=R';
Q=ElectroStaticBalls(XYZ,R, F);
[f, X, Y, p]=SpherePotential(XYZ, Q, R, r0, a, b, Dx, Dy, Nxy);
figure; hold on; grid on; mesh(X,Y,f); 
[Q, D]=ElectroStaticDipoles(XYZ, R, F);
[f, X, Y, p]=SphereDipPotential(XYZ, Q, D, R, r0, a, b, Dx, Dy, Nxy);
figure; hold on; grid on; mesh(X,Y,f); 