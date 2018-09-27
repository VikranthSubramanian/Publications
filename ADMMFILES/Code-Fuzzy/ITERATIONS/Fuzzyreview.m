clc;
clear all
close all
T=0.1; 
alpha=[0.1 0.3 0.5 0.7 0.9];%Threshhold
N=25;    %Number of Grid points
L=600;   %No of samples
Nt=1;    %No of targets
O=[];
O2=[];
O1=[];
tpfac=[];
tpfai=[];
K2=0;
K4=[];
pdc=[];
nvalues=[];
count=0;
Y1=[];
Y2=[];
K6=[];
K7=[];
K8=[];
k11=[];
k12=[];
si=[];
counts=0;
counts1=0;
%% No. of Txs and Rxs
% Mt=4;
% Mr=4;
MtMr=4;
%% Griddefinition
X=[];
x1=[];
Y=[];
y1=[];
for i=1:sqrt(N)
for j=1:sqrt(N)
    x1(j)=j*sqrt(N);  % centroid points
    y1(j)=i*sqrt(N);
end
    X=[X x1];
    Y=[Y y1];
end
X=X';
Y=Y';
TS=[X Y];
X1=[];
%% Tx and Rx Locations and Count
Tx=[7  0;9 6;10 5; 5 4];     %% Tx Positions
Rx=[6  9;9 3;6 8;10 2];     %% Rx Positions
%%  Basis function formulation
counts=0;
counts1=0;
for l=1:length(alpha)
A=ones(1,N);
C1(1,:)=A;
for i=1:(L-1)
C=zeros(size(A));
n=i;  %Shift units
C(n+1:end)=A(1:end-n);
C1(i+1,:)=C;
end
 for i=1:(MtMr)
     si=blkdiag(si,C1);     %% Basis function Matrix
 end
 Y2=[];
 iterations=5;
 for k=1:iterations
  P1=randi(N,1,Nt);
  nvalues=[nvalues;P1];
  s=zeros(1,N*MtMr);
N1=P1;
s(N1)=1;
for i=1:((MtMr)-1)
s(N1+N)=1;
N1=N1+N;
end
for q=-10:5:20
%% Signal Received at the receiver 
z=[];
% q=5;   %%snr    
z=si*s';

 z=awgn(z,q,'measured');
 %% Demodulation at Rx
sdemodulated=pinv(si)*z;
sdemodulated1=reshape(sdemodulated,N,[]);
O=[O;sdemodulated];
for i=1:(MtMr)
 max1(i)=max(sdemodulated1(:,i));
min1(i)=min(sdemodulated1(:,i));
y(:,i)=(sdemodulated1(:,i)-min1(i))/(max1(i)-min1(i));
end
X=[];
Y=[];
[X Y]=fuzzy1(y,alpha(l));
c=X;
count=count+1;
c1(count)=length(c(:,1));
X1=[X1;X];
Y1{count}=Y;
c=0;
Y=0;
end
 K1=find(s);
  for p=1:count
      Y1{p}=reshape(Y1{p},(N*MtMr),[])';
  end
 for j1=1:count
 for i1=1:((MtMr))
     
     X1=Y1{j1};
     Y2=[Y2;Y1];
 if X1(K1(i1))>=T     % Threshold  (Check target position in Rx demodulated signal )
     K2=K2+1; %  (iteration*snrs) 
 end
 X=[];
 K4=[K4;K2];
 K2=0;
 end
 end
 q1=-10:5:20;
SNR1=length(q1);
K4=reshape(K4,MtMr,[]); % Probability of detection
for i2 =1:SNR1
K5(i2)=mean(K4(:,i2));
end
for j=1:count
  for i=1:MtMr
    L1=Y1{j};
     
pfai=find(L1((1+(i-1)*N):(i*N))>0);  % Probability of false alarm indices
pfac=length(pfai);          % No. of false detected points (might include target too)
pf=sdemodulated(pfai);
tpfac=[tpfac;pfac];      %Total count
% tpfai=[tpfai;pfai];        %Total set of indices
O1=[O1;pf];
  end
end
tpfac=reshape(tpfac,MtMr,[]);
tpfac=tpfac-K4;
 for i =1:SNR1
K6(i)=mean(tpfac(:,i));
 end
 K6=K6/(N-1);
 K7=[K7;K6];
 K8=[K8;K5];
 K6=[];
 tpfac=[];
 pfac=[];
 Y2=[Y2;Y1];
%  Y1=[];
 X1=[];
 X=[];
 K1=[];
 K2=0;
 K4=[];
 count=0;
 tpfai=[];
 O1=[O1;O];
 end
 
 for i=1:SNR1
     K9(i)=mean(K7(:,i))/(N-1);
     K10(i)=mean(K8(:,i));
 end
 toc
 counts=counts+1;
 
 subplot(5,2,counts)
 plot(q1,K9);
 subplot(5,2,counts+1)
 plot(q1,K10);
 si=[];
 k11=[k11;K9];
  K9=[];
  k12=[k12;K10];
 K10=[];
 K7=[];
 counts=counts+1;
 K8=[];
 tpfac=[];
 K4=[];
 y=[];
end






