function h=LPN123(f0)
theta=f0(1:2,:);
dtheta=f0(3:4,:)';
gsm=f0(5:7,:);
dgsm=f0(8:10,:);
grn=f0(11:13,:);
dgrn=f0(14:16,:);
lambdasm=f0(17:19,:);
dlambdasm=f0(20:22,:);
lambdarn=f0(23:25,:);
dlambdarn=f0(26:28,:);
p=0.001;
for k=1:5
    theta=theta+(p.*dtheta');
   gsm=gsm+(p.*dgsm);
   grn=grn+(p.*dgrn);
 lambdasm=lambdasm+(p.*dlambdasm);
 lambdarn=lambdarn+(p.*dlambdarn);
end
h=[theta' dtheta dgsm' gsm'  grn' dgrn' lambdasm' dlambdasm' lambdarn' dlambdarn']';
end