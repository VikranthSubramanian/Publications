clc;
clear all
close all
T=0;     %Threshhold
N=25;    %Number of Grid points
L=600;   %No of samples
%% No. of Txs and Rxs
Mt=2;    % No of transmitters
Mr=2;    % No of receivers
%% Initializations
O=[];
O1=[];
tpfac=[];
tpfai=[];
K2=0;
K4=[];
K5=[];
Count12=[];
pdc=[];
si=[];
%% Grid definition
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
%% Tx and Rx Locations and Count
Tx=[7  0;9 6;10 5; 5 4];     %% Tx Positions
Rx=[6  9;9 3;6 8;10 2];     %% Rx Positions
%%  Basis function formulation
A=ones(1,N);
C1(1,:)=A;
for i=1:(L-1)
C=zeros(size(A));
n=i;     %Shift units
C(n+1:end)=A(1:end-n);
C1(i+1,:)=C;
end
 for i=1:(Mt*Mr)
     si=blkdiag(si,C1);     %% Basis function Matrix
 end
%%
Count1=[];
Ans1=[];
Ans=[];
Nvalues=[];
iterations=2;
for k=1:iterations
    P1=randi(N);
        Nvalues=[Nvalues;P1];  %% Accumulation of Random Grid points
   for q=0:2:10                 %%SNR loop
s=zeros(1,N*Mt*Mr);
N1=P1;                       
s(N1)=1;
for i=1:((Mt*Mr)-1)
s(N1+N)=1;
N1=N1+N;
end
%% Signal Received at the receiver 
z=[];
z=si*s';
z=awgn(z,q,'measured');
%% Demodulation at Rx
sdemodulated=pinv(si)*z;
K1=find(s);                   %% Target position index in Transmitted matrix          
for i=1:(Mt*Mr)
if sdemodulated(K1(i))>=T     % Threshold  (Check target position in Rx demodulated signal )
    K2=K2+1;                  %  (iteration*snrs) 
end
K4=[K4;K2];
K5=[K5;K1];
K2=0;
end
% pdc=[pdc;K2];                %Probability detection count 
K2=0;

 for j=1:Mt*Mr
pfai=find(sdemodulated((1+(j-1)*N):(j*N))>0);  % Probability of false alarm indices
pfac=length(pfai);          % No. of false detected points (might include target too)
pf=sdemodulated(pfai);
tpfac=[tpfac;pfac];      %Total count
tpfai=[tpfai;pfai];      %Total set of indices
O1=[O1;pf];
 
 end
O=[O;sdemodulated];
alpha=0.9;
c=reshape(sdemodulated,N,[]);
for i=1:4
    A1{i}=c(:,i);
    end
for i=1:4
 temp=max(A1{i});
 threshhold1=(1-alpha)*temp;    % Deciding the Threshold
 s3=A1{i};
 for j=1:N
    if s3(j)>threshhold1
        Ans=[Ans ;j];
        Count=length(Ans);
    end
  end
    Ans1=[Ans1 ;Ans];
    Ans=[];
    Count1=[Count1;Count]; % In each iteration for each SNR's it gives a Values of No's Grether than threshold.
    Count=0;
end
   end
   S1=reshape(O,N,[]); %% Reshaping for each Sdemodulations
C13=0;
for i=1:(Mt*Mr*(length(0:2:10)))
   C12=Count1(i);
    Ans2=Ans1(C13+1:C12+C13)
    Vals=S1(:,i);
     Vals1=Vals(Ans2);
     ANS3{i}=[Vals1 Ans2]; %% Accumulation of all the Values in pairs like (Mt*Mr*SNR's) for each iteration interations cell array
 C13=C12+C13;
end
ANS4{k}=ANS3;  %% Accumulation of all the Values in pairs like (Mt*Mr*SNR's) in interations cell array
C13=[];
Count12=[Count12;Count1];;
Ans1=[];
Count1=[];  %% Flushing all the values For second iteration
ANS2=[];
S1=[];
A1=[];
c=[];
K4=[];
s=[];
O=[];
z=[];
end











    