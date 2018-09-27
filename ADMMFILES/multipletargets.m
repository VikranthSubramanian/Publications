clc;
clear all
close all
T=0;     %Threshhold
N=25;    %Number of Grid points
L=600;   %No of samples
Nt=2;    %No of targets
%% No. of Txs and Rxs
Mt=4;
Mr=4;
%% Tx and Rx Locations and Count
Tx=[7  0;9 6;10 5; 5 4];     %% Tx Positions
Rx=[6  9;9 3;6 8;10 2];     %% Rx Positions
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
%% Matrix Initialisation 
O=[];
O1=[];
tpfac=[];
tpfai=[];
K2=0;
K4=[];
pdc=[];
si=[];
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
Nvalues=[];
iterations=25;
for i=1:iterations
    P1=randi(N,1,Nt);                      %%% Multiple Targets
        Nvalues=[Nvalues;P1];
   for q=0:3:21                 %%snr loop
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
K1=find(s);  %% Target position index in Transmitted matrix          
for i=1:(Mt*Mr)
if sdemodulated(K1(i))>=T     % Threshold  (Check target position in Rx demodulated signal )
    K2=K2+1; %  (iteration*snrs) 
end
K4=[K4;K2];
K2=0;
end
% pdc=[pdc;K2];                %Probability detection count 
 for j=1:Mt*Mr
pfai=find(sdemodulated((1+(j-1)*N):(j*N))>0);  % Probability of false alarm indices
pfac=length(pfai);          % No. of false detected points (might include target too)
pf=sdemodulated(pfai);
tpfac=[tpfac;pfac];      %Total count
tpfai=[tpfai;pfai];        %Total set of indices
O1=[O1;pf];
 end
 O=[O;sdemodulated];
end
 end
 pdc=K4;
N12=(200);
for i=1:Mt*Mr
    Spd(:,i)=pdc((1+(i-1)*N12):(i*N12));   % Total Probability detection count
end
for i=1:Mt*Mr
    Spfa(:,i)=tpfac((1+(i-1)*N12):(i*N12));  %Total Probability of false alarm
end
Spfa=Spfa-Spd;  % To find the original number of false alarm
for i=1:(Mt*Mr)
    V{i}=reshape(Spd(:,i),length(0:3:21),[]);
    V1{i}=reshape(Spfa(:,i),length(0:3:21),[]);
end
for i=1:(Mt*Mr)
    V{i}=mean(V{i},2);
    V1{i}=mean(V1{i},2)/(N-1);
end
for i=1:Mt*Mr
figure(i);
plot(0:3:21,V{i});    % Plotting Pd probability
figure(i+(Mt*Mr));
plot(0:3:21,V1{i});     % Plotting False alarm probability
end




    
