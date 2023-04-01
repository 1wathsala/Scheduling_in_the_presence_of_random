function [Cost,Position,time] = cmaes(N,C,mu,lmd,p2,p1,cI,cs,cr,VarMin,VarMax,position)
%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA108
% Project Title: Covariance Matrix Adaptation Evolution Strategy (CMA-ES)
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

tic;

%% Problem Settings

%CostFunction=@Ackley;   % Cost Function

nVar=N-1;                % Number of Unknown (Decision) Variables

VarSize=[1 nVar];       % Decision Variables Matrix Size


%% CMA-ES Settings

% Maximum Number of Iterations
MaxIt=5;

% Population Size (and Number of Offsprings)
lambda=(4+round(3*log(nVar)))*10;

% Number of Parents
mut=round(lambda/2);

% Parent Weights
w=log(mut+0.5)-log(1:mut);
w=w/sum(w);

% Number of Effective Solutions
mu_eff=1/sum(w.^2);

% Step Size Control Parameters (c_sigma and d_sigma);
sigma0=0.3*(VarMax-VarMin);
cs1=(mu_eff+2)/(nVar+mu_eff+5);
ds=1+cs1+2*max(sqrt((mu_eff-1)/(nVar+1))-1,0);
ENN=sqrt(nVar)*(1-1/(4*nVar)+1/(21*nVar^2));

% Covariance Update Parameters
cc=(4+mu_eff/nVar)/(4+nVar+2*mu_eff/nVar);
c1=2/((nVar+1.3)^2+mu_eff);
alpha_mu=2;
cmu=min(1-c1,alpha_mu*(mu_eff-2+1/mu_eff)/((nVar+2)^2+alpha_mu*mu_eff/2));
hth=(1.4+2/(nVar+1))*ENN;

%% Initialization

ps=cell(MaxIt,1);
pc=cell(MaxIt,1);
C1=cell(MaxIt,1);
sigma=cell(MaxIt,1);

ps{1}=zeros(VarSize);
pc{1}=zeros(VarSize);
C1{1}=eye(nVar);
sigma{1}=sigma0;

empty_individual.Position=[];
empty_individual.Step=[];
empty_individual.Cost=[];

M=repmat(empty_individual,MaxIt,1);
%M(1).Position=unifrnd(VarMin,VarMax,VarSize);
M(1).Position = position;
M(1).Step=zeros(VarSize);
%M(1).Cost=CostFunction(M(1).Position);
M(1).Cost = cost_fun_eval(M(1).Position,N,C,mu,lmd,p2,p1,cI,cs,cr);

BestSol=M(1)

BestCost=zeros(MaxIt,1);

%% CMA-ES Main Loop

for g=1:MaxIt
    
    % Generate Samples
    pop=repmat(empty_individual,lambda,1);
    parfor i=1:lambda
        pop(i).Step=mvnrnd(zeros(VarSize),C1{g});
        pop(i).Position=M(g).Position+sigma{g}*pop(i).Step;
        if pop(i).Position >= VarMin & pop(i).Position <= VarMax
            pop(i).Cost=cost_fun_eval(pop(i).Position,N,C,mu,lmd,p2,p1,cI,cs,cr);
        
            % Update Best Solution Ever Found
        else
            pop(i) = BestSol
        end
    end
    
    % Sort Population
    Costs=[pop.Cost];
    [Costs, SortOrder]=sort(Costs);
    pop=pop(SortOrder);
  
    if pop(1).Cost<BestSol.Cost
        BestSol=pop(1);
    end
    % Save Results
    BestCost(g)=BestSol.Cost;
    
    % Display Results
    disp(['Iteration ' num2str(g) ': Best Cost = ' num2str(BestCost(g))]);
    BestSol
    
    % Exit At Last Iteration
    if g==MaxIt
        break;
    end
        
    % Update Mean
    M(g+1).Step=0;
    for j=1:mut
        M(g+1).Step=M(g+1).Step+w(j)*pop(j).Step;
    end
    M(g+1).Position=M(g).Position+sigma{g}*M(g+1).Step;
    if M(g+1).Position >= VarMin & M(g+1).Position <= VarMax
        M(g+1).Cost=cost_fun_eval(M(g+1).Position,N,C,mu,lmd,p2,p1,cI,cs,cr);
        if M(g+1).Cost<BestSol.Cost
            BestSol=M(g+1);
        end
    else
        M(g+1) = BestSol;
    end
    
    % Update Step Size
    ps{g+1}=(1-cs1)*ps{g}+sqrt(cs1*(2-cs1)*mu_eff)*M(g+1).Step/chol(C1{g})';
    sigma{g+1}=sigma{g}*exp(cs1/ds*(norm(ps{g+1})/ENN-1))^0.3;
    
    % Update Covariance Matrix
    if norm(ps{g+1})/sqrt(1-(1-cs1)^(2*(g+1)))<hth
        hs=1;
    else
        hs=0;
    end
    delta=(1-hs)*cc*(2-cc);
    pc{g+1}=(1-cc)*pc{g}+hs*sqrt(cc*(2-cc)*mu_eff)*M(g+1).Step;
    C1{g+1}=(1-c1-cmu)*C1{g}+c1*(pc{g+1}'*pc{g+1}+delta*C1{g});
    for j=1:mut
        C1{g+1}=C1{g+1}+cmu*w(j)*pop(j).Step'*pop(j).Step;
    end
    
    % If Covariance Matrix is not Positive Defenite or Near Singular
    [V, E]=eig(C1{g+1});
    if any(diag(E)<0)
        E=max(E,0);
        C1{g+1}=V*E/V;
    end
      
end

Cost = BestSol.Cost
Position = BestSol.Position

%% Display Results

%figure;
% plot(BestCost, 'LineWidth', 2);
%semilogy(BestCost, 'LineWidth', 2);
%xlabel('Iteration');
%ylabel('Best Cost');
%grid on;
time = toc;
end

