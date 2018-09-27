clc;
clear;
close all;
% % % m=1;
% n=100;
wmax=0.9;
wmin=0.4;
c1=2;
c2=2;
% maxite=10;
N=25;    %Number of Grid points
L=625;   %No of samples
%% No. of Txs and Rxs
Mt=4;
Mr=4;
si=[];
mse=[];
%%  Basis function formulation
% A=ones(1,N);
% C1(1,:)=A;
% for i=1:(L-1)
% C=zeros(size(A));
% n=i;     %Shift units
% C(n+1:end)=A(1:end-n);
% C1(i+1,:)=C;
% end
%  for i=1:(4)
%      si=blkdiag(si,C1);     %% Basis function Matrix
%  end
A=ones(1,125);
B=zeros(1,(L-125));
C=[A B]';
C12=[];
C1=[];
for v=1:N
C12=circshift(C,(v-1));
C1=[C1 C12];
end
 for i=1:(Mt*Mr)
     si=blkdiag(si,C1);     %% Basis function Matrix
 end
  P1=randi(N);
  s=zeros(1,N*Mt*Mr);
N1=P1;
s(N1)=1;
for i=1:((Mt*Mr)-1)
s(N1+N)=1;
N1=N1+N;
end
BestSol1=[];
BestSol2=[];
BestSol3=[];
countz=0;
for q=-10:5:20
z=si*s';
 z=awgn(z,q,'measured');
%% Demodulation at Rx
sdemodulated=pinv(si)*z;
% sdemodulated=sdemodulated';
% for i=1:(Mt*Mr)
 max1=max(sdemodulated);
min1=min(sdemodulated);
y=(sdemodulated-min1)/(max1-min1);

%% Problem Definition

CostFunction=@(x,si,z) Sphere(x,si,z);        % Cost Function

nVar=(N*Mt*Mr);            % Number of Decision Variables

VarSize=[1 nVar];   % Size of Decision Variables Matrix

VarMin=0;         % Lower Bound of Variables
VarMax= 1;         % Upper Bound of Variables


%% PSO Parameters

MaxIt=1250;      % Maximum Number of Iterations

nPop=100;        % Population Size (Swarm Size)

% PSO Parameters
w=1;            % Inertia Weight
wdamp=0.99;     % Inertia Weight Damping Ratio
c1=1.5;         % Personal Learning Coefficient
c2=2.0;         % Global Learning Coefficient

% If you would like to use Constriction Coefficients for PSO,
% uncomment the following block and comment the above set of parameters.

% % Constriction Coefficients
% phi1=2.05;
% phi2=2.05;
% phi=phi1+phi2;
% chi=2/(phi-2+sqrt(phi^2-4*phi));
% w=chi;          % Inertia Weight
% wdamp=1;        % Inertia Weight Damping Ratio
% c1=chi*phi1;    % Personal Learning Coefficient
% c2=chi*phi2;    % Global Learning Coefficient

% Velocity Limits
VelMax=0.1*(VarMax-VarMin);
VelMin=-VelMax;

%% Initialization

empty_particle.Position=[];
empty_particle.Cost=[];
empty_particle.Velocity=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];

particle=repmat(empty_particle,nPop,1);

GlobalBest.Cost=inf;
% P2=randi(N);
for i=1:nPop
    Y{i}=rand(1,N*Mt*Mr);
    N2=P1;
    Y1=Y{i};
    Y1(N2)=1;
for j=1:((Mt*Mr)-1)
Y1(N2+N)=1;
N2=N2+N;
end
Y2{i}=Y1;
end

for i=1:nPop
    
    % Initialize Position
   particle(i).Position=Y2{i};
    
    % Initialize Velocity
    particle(i).Velocity=zeros(VarSize);
    
    % Evaluation
    particle(i).Cost=CostFunction(particle(i).Position,si,z);
    
    % Update Personal Best
    particle(i).Best.Position=particle(i).Position;
    particle(i).Best.Cost=particle(i).Cost;
    
    % Update Global Best
    if particle(i).Best.Cost<GlobalBest.Cost
        
        GlobalBest=particle(i).Best;
        
    end
    
end

BestCost=zeros(MaxIt,1);

%% PSO Main Loop

for it=1:MaxIt
    
    for i=1:nPop
        
        % Update Velocity
        particle(i).Velocity = w*particle(i).Velocity ...
            +(c1*rand(VarSize)).*(particle(i).Best.Position-particle(i).Position) ...
            +(c2*rand(VarSize)).*(GlobalBest.Position-particle(i).Position);
        
        % Apply Velocity Limitsparticle(i).Position
        particle(i).Velocity = max(particle(i).Velocity,VelMin);
        particle(i).Velocity = min(particle(i).Velocity,VelMax);
        
        % Update Position
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        
        % Velocity Mirror Effect
        IsOutside=(particle(i).Position<VarMin | particle(i).Position>VarMax);
        particle(i).Velocity(IsOutside)=-particle(i).Velocity(IsOutside);
        
        % Apply Position Limits
        particle(i).Position = max(particle(i).Position,VarMin);
        particle(i).Position = min(particle(i).Position,VarMax);
        
        % Evaluation
        particle(i).Cost = CostFunction(particle(i).Position,si,z);
        
        % Update Personal Best
        if particle(i).Cost<particle(i).Best.Cost
            
            particle(i).Best.Position=particle(i).Position;
            particle(i).Best.Cost=particle(i).Cost;
            
            % Update Global Best
            if particle(i).Best.Cost<GlobalBest.Cost
                
                GlobalBest=particle(i).Best;
                
            end
            
        end
      
    end
    
BestCost(it)=GlobalBest.Cost;
    
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    
    w=w*wdamp;
    
end
countz=countz+1;
BestSol = GlobalBest;
% % BestSol2=BestSol.Position;
 mse(countz)=norm(BestSol.Position-s)/(N*Mt*Mr);
% mse(countz)=norm(BestSol.Cost)/100;

BestSol1=[BestSol1 BestSol];
% BestSol2=[];
end

%% Results

% figure;
% %plot(BestCost,'LineWidth',2);
% plot(BestCost,'LineWidth',2);
% xlabel('Iteration');
% ylabel('Best Cost');
% grid on;
plot(-10:5:20,mse);