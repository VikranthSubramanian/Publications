function [u,c23] =fuzzy1(yt,MtMr1)
y2=yt;
alpha=0.1;
c=y2;
c1=[];
N=25;
counts=0;
counts1=0;
MtMr=MtMr1;
Ans=[];
Ans1=[];
Count=0;
Count1=[];
for i=1:(MtMr)
    A1{i}=c(:,i);
    end
for i=1:(MtMr)
 temp=max(A1{i});
 threshhold1=(1-alpha)*temp;    % Deciding the Threshold
 s3=A1{i};
 for j=1:N
    if s3(j)<threshhold1
        s3(j)=0;
    end
  end
    Ans1=[Ans1 ;s3];
    Ans=[];
    Count1=[Count1;Count]; % In each iteration for each SNR's it gives a Values of No's Grether than threshold.
    Count=0;
end
Ans1=reshape(Ans1,N,[]);
   c23=Ans1;
   
for i=1:(MtMr)
    [c(:,i),ind1(:,i)]=sort(c(:,i));
    c1(:,i)=c(N-1:N,i);
     ind(:,i)=ind1(N-1:N,i);
end
count=0;
u(:,1)=unique(ind);
ind=reshape(ind,1,[]);
for w=1:length(u)
    for e=1:length(ind)
        if u(w)==ind(e)
            count=count+1;
        end
    end
    u(w,2)=count;
    count=0;
end
s3=0;
Ans1=[];
end