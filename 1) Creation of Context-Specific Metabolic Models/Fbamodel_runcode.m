%load('geni.mat')
%load('reaction_expression.mat')
%load('Microarraydata829.mat') 
%load('rCAI_norm829.mat')

addpath(genpath('.\cobra'));
initCobraToolbox

M=2;
V = numel(geni);

Objective_1=0;
Objective_2=0;


for i=1:size(Microarraydata834,2)
    %x=ones(834,1).*rCAI_norm834; 1.47
    x=Microarraydata834(:,i).*rCAI_norm834;
    tmp = evaluate_objective(x,2,numel(geni),fbamodel,geni,reaction_expression);
    Objective_1(i)=-tmp(1);
    Objective_2(i)=tmp(2);
end


