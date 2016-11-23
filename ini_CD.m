function [D1,D2]= ini_CD(Nz,dz)

% Differentiation matrices of the compact scheme
R= [1 1/3 zeros(1,Nz-2) 1/3];
C= [1 1/3  zeros(1,Nz-2) 1/3];
A= toeplitz(C,R);
A(1,:)= [1 2 zeros(1,Nz-1)];
A(2,:)= [1/4 1 1/4 zeros(1,Nz-2)];
A(end-1,:)= [zeros(1,Nz-2) 1/4 1 1/4];
A(end,:)= [zeros(1,Nz-1) 2 1];

R= [0 7/9 1/36 zeros(1,Nz-4) -1/36 -7/9]/dz;
C= [0 -7/9 -1/36  zeros(1,Nz-4) 1/36 7/9]/dz;
B= toeplitz(C,R);
B(1,:)= [-5/2 2 1/2 zeros(1,Nz-2)]/dz;
B(2,:)= [-3/4 0 3/4 zeros(1,Nz-2)]/dz;
B(end-1,:)= [zeros(1,Nz-2) -3/4 0 3/4]/dz;
B(end,:)= [zeros(1,Nz-2) -1/2 -2 5/2]/dz;


D1= A\B;
D2= D1*D1;

A(end,:)= 0;A(end,end)= 1;
B(end,:)= 0;
D1= A\B;