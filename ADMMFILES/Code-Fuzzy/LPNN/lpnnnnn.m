%   function f=lpnnnnn(t0,f0)
% C0=0.00000001; 
C0=5;
p=0.001; % Stepsize
sigma =10; % Noise Variance
%No. of Txs and Rxs
M=3;
N=3;
c=3;
%Tx and Rx positions
Tx=[-1000 -1300 ;    
     500 2000 ;
     2500 0];
Rx=[1500 -1800 ;
   2100 1500; 
   -1200 1000];
%%Variables=  theta,gsm,grn,lambdasm,lambdarn
%% initialisation
b=0.0000005;
a=0.0002;
t=(b-a).*rand(M,N)+a;
theta=rand(2,1)';
gsm=12*rand(1,M)';
lambdasm=rand(1,M)';
grn=200*rand(1,N)';
lambdarn=rand(1,N)';
% theta=1200*f0(1:2,:)
% gsm=12*f0(3:5,:)
% grn=1000*f0(6:8,:)
% lambdasm=f0(9:11,:)
% lambdarn=f0(12:14,:)
% dtheta=zeros(1,2);
theta1=[0 200];
dgsm=[200;200;200];
dgrn=[200;200;200];
dlambdasm=zeros(M,1);
dlambdarn=zeros(N,1);
%% dtheta
dtheta1=[0 0];
dtheta2=[0 0];
% Problem in dtheta first sought it out. Then only all the values will be
% corrected
for i=1:M
   dtheta11(i,:)=2*(lambdasm(i)+(C0*((gsm(i)^2)-((theta(1,:)-Tx(i,:))*(theta(1,:)-Tx(i,:))'))))*(theta(1,:)-Tx(i,:)) ;
   dtheta1(1,:)=dtheta1(1,:)+dtheta11(i,:);
end
for j=1:N
   dtheta21(j,:)=2*(lambdarn(j)+(C0*((grn(j)^2)-((theta(1,:)-Rx(j,:))*(theta(1,:)-Rx(j,:))'))))*(theta(1,:)-Rx(j,:)) ;
   dtheta2(1,:)=dtheta2(1,:)+dtheta21(j,:);
end
dtheta=dtheta1+dtheta2
% bistatic measurement matrix
for i=1:M
   for j=1:N
       q(i,j)=((theta1(1,:)-Tx(i,:))*(theta1(1,:)-Tx(i,:))')+((theta1(1,:)-Rx(j,:))*(theta1(1,:)-Rx(j,:))');
   end
end
q=t*c;

 r = awgn(q,20,'measured');
%% dgsm
c1=0;
for i=1:M
    for j=1:N
   c1=c1+(r(i,j)-grn(j));
   end
r2(i)=c1
end
for i=1:M
dgsm(i)=(2/(sigma))*(r2(i)-gsm(i))-(2*gsm(i)*lambdasm(i))-(2*C0*gsm(i)*((gsm(i)^2)-((theta(1,:)-Tx(i,:))*(theta(1,:)-Tx(i,:))')))
end

%% dgrn
c2=0;
for j=1:N
    for i=1:M
   c2=c2+(r(j,i)-gsm(i));
   end
r3(j)=c2;
end
for j=1:N
dgrn(j)=((2/(sigma))*(r3(j)-grn(j))-(2*grn(j)*lambdarn(j))-(2*C0*grn(j)*((grn(j)^2)-((theta(1,:)-Rx(j,:))*(theta(1,:)-Rx(j,:))'))));
end
%% dlambdasm
for i=1:M
    dlambdasm(i)=((gsm(i)^2)-((theta(1,:)-Tx(i,:))*(theta(1,:)-Tx(i,:))'));
end
%% dlambdarn
for j=1:N
    dlambdarn(j)=((grn(j)^2)-((theta(1,:)-Rx(j,:))*(theta(1,:)-Rx(j,:))'));
 end
% f=[theta' gsm' grn' lambdasm' lambdarn']'
f=[theta dtheta dgsm' gsm'  grn' dgrn' lambdasm' dlambdasm' lambdarn' dlambdarn']'
g=LPN123(f);
 
%   end