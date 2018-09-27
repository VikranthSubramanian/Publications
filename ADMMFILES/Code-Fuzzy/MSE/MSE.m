
T=0;     %Threshhold
N=25;    %Number of Grid points
L=625;   %No of samples
%% No. of Txs and Rxs
Mt=4;
Mr=4;
%% 
O=[];
si=[];
A=ones(1,125);
B=zeros(1,(L-125));
C=[A B]';
C12=[];
mse1=[];
mse2=[];
C1=[];
Nvalues=[];
count=0;
iterations=1;
for v=1:N
C12=circshift(C,(v-1));
C1=[C1 C12];
end

 for i=1:(Mt*Mr)
     si=blkdiag(si,C1);     %% Basis function Matrix
 end
 for i=1:iterations
    P1=randi(N);
        Nvalues=[Nvalues;P1];
   for q=-10:3:10                 %%snr loop
s=zeros(1,N*Mt*Mr);
N1=P1;
s(N1)=1;
for i=1:((Mt*Mr)-1)
s(N1+N)=1;
N1=N1+N;
end
z=si*s';
 z=awgn(z,q,'measured');
%% Demodulation at Rx
sdemodulated=pinv(si)*z;
count=count+1;
mse(count)=norm(sdemodulated-s')/(N*Mt*Mr);
   end
   mse1=[mse1 mse];
   mse=[];
   count=0;
 end
mse1=reshape(mse1,iterations,[]);
for i=1:7
    mse2(i)=mean(mse1(:,i));
end
plot(-10:3:10,mse2);

