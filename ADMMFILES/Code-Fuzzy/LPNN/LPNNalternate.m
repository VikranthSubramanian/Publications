clc;
clear all;
close all;
%% Tx and Rx positions
M=3;
N=3;
sm=[-1000 500  2500
    -1300 2000 0];
rn=[1500  2100 -1200 
   -1800 1500  1000];
   
%% Initialisation
C0=5;
p=0.001;     % Stepsize
sigma=10;    %Noise variance
iterations=50;
I=eye(2);
I1=zeros(M+N,M+N);
oneM=ones(M,1);
oneN=ones(N,1);
theta0=[200 
    200];
gsm0=rand(M,1); 
grn0=rand(N,1);
lamdasm0=rand(M,1); 
lamdarn0=rand(N,1);
l1=0;l3=0;
l2=zeros(2,1);
l4=zeros(2,1);
theta=[];
gsm=[]; 
grn=[];
lamdasm=[]; 
lamdarn=[];
%% Equilibrium points 
% thetae=zeros(2,1);
% gsme=zeros(M,1);
% grne=zeros(N,1);
% lamdasme=zeros(M,1);
% lamdarne=zeros(N,1);
thetae=[0
    200];
gsme=zeros(M,1);
grne=zeros(N,1);
lamdasme=zeros(M,1);
lamdarne=zeros(N,1);
%% 
for i=1:M
    d1(:,i)=(thetae(:,1)-sm(:,i));
end
for j=1:N
    d2(:,j)=(thetae(:,1)-rn(:,j));
end
for i=1:M
    l1=l1+lamdasme(i);
end
for j=1:N
    l3=l3+lamdarne(j);
end
a11=((-2)*(l1+l3)*I)+(4*C0*0.00001*((d1*d1')+(d2*d2')));
for i=1:M
a12(:,i)=(-4*C0)*(gsme(i)*(thetae(:,1)-sm(:,i)));
end
for j=1:N
a13(:,j)=(-4*C0)*(grne(j)*(thetae(:,1)-rn(:,j)));
end
a22=(2*diag(lamdasme))+ (4*C0*(gsme'*gsme)*eye(M))+(((2*N)/sigma)*eye(M));
a23=((2/sigma)*(oneM*oneN'));
a33=(2*diag(lamdarne))+ (4*C0*(grne'*grne)*eye(N))+(((2*M)/sigma)*eye(N));
a21=-a12';
a31=-a13';
a32=-a23';
for i=1:M
b11(:,i)=(-2)*(thetae(:,1)-sm(:,i));
end
for j=1:N
b12(:,j)=(-2)*(thetae(:,1)-rn(:,j));
end
b21=2*diag(gsme);
b22=2*diag(grne);
b31=zeros(M,M);
b32=zeros(N,N);
a1=[a11 a12 a13];
a2=[a21 a22 a23];
a3=[a31 a32 a33];
a=[a1
    a2
    a3];
b1=[b11 b12];
b2=[b21 b22];
b3=[b31 b32];
b=[b1
    b2
    b3];
G=-[a b
    -b' I1];


%% Iterations
for k=1:iterations
    dummy=[theta0-thetae
    gsm0-gsme
    grn0-grne
    lamdasm0-lamdasme
    lamdarn0-lamdarne];
f=-G*dummy;
dtheta=f(1:2,1);
dgsm=f(3:5,1);
dgrn=f(6:8,1);
dlamdasm=f(9:11,1);
dlamdarn=f(12:14,1);  
theta1=theta0-(p*dtheta);
        gsm1=gsm0+(p*dgsm);    
        grn1=grn0+(p*dgrn);
            lamdasm1=lamdasm0+(p*dlamdasm);
                lamdarn1=lamdarn0+(p*dlamdarn);
                theta=[theta theta1];
                gsm=[gsm gsm1];
                grn=[grn grn1];
                lamdasm=[lamdasm lamdasm1];
                lamdarn=[lamdarn lamdarn1];
                theta0=theta1;
                gsm0=gsm1;
                grn0=grn1;
                lamdasm0=lamdasm1;
                lamdarn0=lamdarn1;
end