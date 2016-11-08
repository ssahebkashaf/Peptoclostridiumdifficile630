load('geni.mat')
load('reaction_expression.mat')
%load('Microarraydata834.mat')
load('rCAI_norm834.mat')



M=2;
V = numel(geni);
%x = ones(1,numel(geni)); %load x which I obtained from microarray analysis
%x=x'

Objective_1=0;
Objective_2=0;

%x = ones(numel(geni),1);%.*rCAI_norm834*10;%Microarraydata834(:,i)
x = x1.*rCAI_norm834
obj_unperturbed = evaluate_objective(x,2,numel(geni),fbamodel,geni,reaction_expression);

% 
% b=0
% j=1
% for i=1:length(val_pathway)  %make val_pathway from SS_gene
% if val_pathway(i)==1
% b(j)=i
% j=j+1
% end
% end

mca = zeros(numel(geni),1);
mca2 = zeros(numel(geni),1);

for i=1:numel(geni)
    y = x;
    y(i) = y(i) - 0.005*10;
    %y = y.*rCAI_norm834;
    
    tmp = evaluate_objective(y,2,numel(geni),fbamodel,geni,reaction_expression);
    tmp1(i)=tmp(1);
    tmp2(i)=tmp(2);
end
